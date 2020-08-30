#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates fasta file out of provided peaks'''

import argparse
import os
import sys
from collections import defaultdict, Counter


from Bio import SeqIO

import numpy as np;
from pybedtools import BedTool, Interval
from afbio.pybedtools_af import peak2seq



parser = argparse.ArgumentParser(description='Creates fasta file out of provided peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the peaks/regions, gff/bed format");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file")
parser.add_argument('--flank', nargs = '?', required=True, type = int, help = "Flank length of sequence")
parser.add_argument('--mincov', nargs = '?', default=0.0, type = float, help = "Minimum required peak coverage")
args = parser.parse_args();

genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
for interval in BedTool(args.path):
    covscore = np.mean([float(x) for x in interval.attrs['topcoverage'].split(',')])
    if(covscore >= args.mincov):
        seq = peak2seq(interval, genome, args.flank)
        sys.stdout.write( ">%s:%s:%d:%d\n%s\n" % (interval.chrom, interval.strand, interval.start, interval.end, seq))
    

            
            
