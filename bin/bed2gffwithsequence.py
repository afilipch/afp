#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Adds sequence annotation to the provided bed file'''

import argparse
import os
import sys
from pybedtools import BedTool
from Bio import SeqIO

from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Adds sequence annotation to the provided bed file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the bed file");
parser.add_argument('--genome', nargs = '?', required = True, type = str, help = "Path to the genome, fasta format");
parser.add_argument('--gff', nargs = '?', default = False, const = True, type = bool, help = "If set, output is in gff format");
args = parser.parse_args();

genome = next(SeqIO.parse(args.genome, 'fasta')).seq

              
for interval in BedTool(args.path):
    seq = str(genome[interval.start: interval.end].upper())
    if(args.gff):
        sys.stdout.write(str(construct_gff_interval(interval.chrom, interval.start, interval.end, 'un', score=interval.score, strand=interval.strand, source='un', frame='.', attrs=[('Name', interval.name), ('seq', seq)] )))
    else:
        print("\t".join([str(x) for x in list(interval) + [seq]]));
