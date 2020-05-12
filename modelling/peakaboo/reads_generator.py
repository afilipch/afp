#! /usr/local/anaconda3/bin/python
'''Generates artificial chap-seq reads'''

import argparse
import os
import sys
import numpy as np;
import random

from Bio import SeqIO
from pybedtools import BedTool

from afbio.generators import generator_doublesam
from afbio.numerictools import select_by_probability

parser = argparse.ArgumentParser(description='Generates artificial chap-seq reads');

parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to genome, fasta file");
parser.add_argument('--peaks', nargs = '?', required=True, type = str, help = "Path to the real peaks, gff format");

parser.add_argument('--numreads', nargs = '?', default=400000, type = int, help = "Number of reads to generate");
parser.add_argument('--numsamples', nargs = '?', default=50, type = int, help = "Number of samples to generate");
parser.add_argument('--lengths', nargs = '+', default=[300], type = int, help = "Length of the generated reads");
parser.add_argument('--bs_width', nargs = '?', default=20, type = int, help = "Size of bindig sites");
parser.add_argument('--fraction_peaks', nargs = '?', default=0.2, type = float, help = "Fraction of reads which belong to peaks");
parser.add_argument('--errors', nargs = '+', default=[0], type = int, help = "Sequence error rate in integer percents");

parser.add_argument('--mincov', nargs = '?', default=3, type = float, help = "Minimum strength of the peaks");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();



NUCLS = 'ACTG'
FACTOR = len(NUCLS)/((len(NUCLS) - 1)*100)





def mutate_sequence(seq, mmrate):
    res = []
    for s in seq:
        if(random.random() < mmrate):
            res.append(random.choice(NUCLS))
        else:
            res.append(s)
    return ''.join(res)
    

def get_nonpeak_read(genome, chrnames, chrsizes, length, mmrate):
    chrom = select_by_probability(chrnames, chrsizes);
    seqrecord = genome[chrom];
    ul = len(seqrecord) - length
    start = random.randint(0, ul);
    local = seqrecord[start:start+length]
    if(random.randint(0,1)):
        seq = str(local.seq.upper())
        strand = '+'
    else:
        seq = str(local.seq.reverse_complement().upper())
        strand = '-'
    if(mmrate):
        seq = mutate_sequence(seq, mmrate)
    return seq, peak.chrom, start, start+length, strand

    
    
def get_peak_read(peak, length, genome, mmrate):
    start = random.randint(peak.stop - length, peak.start)
    local = genome[peak.chrom][start:start+length];
    if(peak.strand == '+'):
        seq = str(local.seq.upper())
        strand = '+'
    else:
        seq = str(local.seq.reverse_complement().upper())
        strand = '-'
    if(mmrate):
        seq = mutate_sequence(seq, mmrate)
    return seq, peak.chrom, start, start+length, strand
    



def convert(a):
    header = "|".join([ str(x) for x in a[1:]])
    return ">%s\n%s" % (header, a[0])


genome = SeqIO.to_dict(SeqIO.parse(args.genome, 'fasta'));
chrnames = list(genome.keys())
chrsizes = [len(genome[x]) for x in chrnames] 

peaks = [x for x in BedTool(args.peaks) if float(x.score) > args.mincov]
for peak in peaks:
    peak.start = int(peak.name) - int(args.bs_width/2)
    peak.stop = int(peak.name) + int(args.bs_width/2)
peaks_probabilites = [float(peak.score) for x in peaks]




for length in args.lengths:
    for error in args.errors:
        for sc in range(args.numsamples):
            with open(os.path.join(args.outdir, "sample_%d_%d_%d.fa" % (length, error, sc+1)), 'w') as f:
                mmrate = FACTOR*error
                for _ in range(args.numreads):
                    if(random.random()<=args.fraction_peaks):
                        selected_peak = select_by_probability(peaks, peaks_probabilites);
                        f.write("%s\n" %  convert(get_peak_read(selected_peak, length, genome, mmrate)) ) 
                    else:
                        f.write("%s\n" % convert(get_nonpeak_read(genome, chrnames, chrsizes, length, mmrate)) )





    
    
