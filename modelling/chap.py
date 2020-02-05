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

parser.add_argument('--numreads', nargs = '?', default=1000000, type = int, help = "Number of reads to generate");
parser.add_argument('--length', nargs = '?', default=300, type = int, help = "Length of the generated reads");
parser.add_argument('--bs_width', nargs = '?', default=20, type = int, help = "Size of bindig sites");
parser.add_argument('--fraction_peaks', nargs = '?', default=0.2, type = int, help = "Fraction of reads which belong to peaks");

parser.add_argument('--mincov', nargs = '?', default=3, type = float, help = "Minimum strength of the peaks");
#parser.add_argument('--num_peaks', nargs = '?', default=300, type = int, help = "Number of peaks to generate");
args = parser.parse_args();

def get_nonpeak_read(genome, chrnames, chrsizes, length):
    chrom = select_by_probability(chrnames, chrsizes);
    seqrecord = genome[chrom];
    ul = len(seqrecord) - length
    start = random.randint(0, ul);
    local = seqrecord[start:start+length]
    if(random.randint(0,1)):
        return str(local.seq.upper()), chrom, start, start+length, '+'
    else:
        return str(local.seq.reverse_complement().upper()), chrom, start, start+length, '-'
    
    
def get_peak_read(peak, length, genome):
    start = random.randint(peak.stop - length, peak.start)
    local = genome[peak.chrom][start:start+length];
    if(peak.strand == '+'):
        return str(local.seq.upper()), peak.chrom, start, start+length, '+'
    else:
        return str(local.seq.reverse_complement().upper()), peak.chrom, start, start+length, '-'
    


def convert(a):
    header = "|".join([ str(x) for x in a[1:]])
    return "@%s\n%s" % (header, a[0])


genome = SeqIO.to_dict(SeqIO.parse(args.genome, 'fasta'));
chrnames = list(genome.keys())
chrsizes = [len(genome[x]) for x in chrnames] 


peaks = [x for x in BedTool(args.peaks) if float(x.attrs['topcoverage']) > args.mincov]
for peak in peaks:
    peak.start = int(peak.name) - int(args.bs_width/2)
    peak.stop = int(peak.name) + int(args.bs_width/2)
    
    
peaks_probabilites = [float(x.attrs['topcoverage']) for x in peaks]

for _ in range(args.numreads):
    if(random.random()<=args.fraction_peaks):
        selected_peak = select_by_probability(peaks, peaks_probabilites);
        print(convert(get_peak_read(selected_peak, args.length, genome))) 
    else:
        print(convert(get_nonpeak_read(genome, chrnames, chrsizes, args.length)))





    
    
