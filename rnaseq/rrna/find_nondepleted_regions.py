#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Detects nondepleted rRNA regions beased on the genomic coverage'''

import argparse
import os
import sys
from collections import defaultdict

from pybedtools import BedTool
import numpy as np
from Bio import SeqIO

from afbio.generators import get_only_files
from afbio.sequencetools import coverage2dict

parser = argparse.ArgumentParser(description='Detects nondepleted rRNA regions beased on the genomic coverage');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the coverage folder");
parser.add_argument('--rrna', nargs = '?', required=True, type = str, help = "Path to the rRNA, gff file");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
parser.add_argument('--minfraction', nargs = '?', default=1, type = float, help = "Minimal required fraction/multiplier (of the mean rRNA coverage) for a particular position to be counted as nondepleted");
parser.add_argument('--minlength', nargs = '?', default=20, type = float, help = "Minimal required length of non-depleted regions");
args = parser.parse_args();

strand_conv = {'plus': '+', 'minus': '-'} 

genome = SeqIO.to_dict(SeqIO.parse(args.genome,'fasta'))

#rrna = BedTool(args.rrna);
files = get_only_files(args.path)
strand2coverage = {}
for f in files[:]:
    strand = strand_conv[f.split(".")[-2]]
    for chrom, cov in coverage2dict(f, cpos=2).items():
        if( (chrom, strand) in strand2coverage):
            strand2coverage[(chrom, strand)] += cov
        else:
            strand2coverage[(chrom, strand)] = cov
            

total_sum = sum([sum(x) for x in strand2coverage.values()])
#print(total_sum)

rrna2coverage =[];
for interval in BedTool(args.rrna):
    local_coverage = strand2coverage[(interval.chrom, interval.strand)][interval.start:interval.stop]
    rrna2coverage.append(( interval.chrom, interval.start, interval.stop, interval.attrs['Parent'], interval.strand, interval.attrs['product'], local_coverage ))
    
rrna_sum = sum([ sum(x[-1]) for x in rrna2coverage ])
rrna_length = sum([ x[2]-x[1] for x in rrna2coverage ])
covlimit = args.minfraction*(rrna_sum/rrna_length)



def print_out(chrom, r_start, r_stop, cur_start, pos, strand, name, type_, cov, genome):
    start = r_start + cur_start
    stop = r_start + pos
    seq = genome[chrom][start: stop].seq
    if(strand == '-'):
        seq = seq.reverse_complement()
    print( "\t".join(( chrom, str(start), str(stop), strand, name, type_, "%1.1f" % np.mean(cov[cur_start: pos]), str(seq) )) )
    

print( "\t".join(( "chromosome", "start", "stop", "strand", "rRNA_ID", "type", "mean_coverage", 'sequence' )) )
for chrom, r_start, r_stop, name, strand, type_, cov in rrna2coverage:
    in_region = cov[0] >= covlimit
    cur_start = 0;
    for pos, c in enumerate(cov[1:], start=1):
        if(c >= covlimit):
            if(not in_region):
                cur_start = pos;
                in_region = True;
        else:
            if(in_region):
                in_region = False;
                if(pos - cur_start >= args.minlength):
                    print_out(chrom, r_start, r_stop, cur_start, pos, strand, name, type_, cov, genome)
    else:
        if(in_region):
            pos = len(cov)
            if(pos - cur_start >= args.minlength):
                print_out(chrom, r_start, r_stop, cur_start, pos, strand, name, type_, cov, genome)
        
                
    


#plus_files = [x for x in files if x.split(".")[-2] == 'plus']
#minus_files = [x for x in files if x.split(".")[-2] == 'minus']
#print(plus_files)



        

