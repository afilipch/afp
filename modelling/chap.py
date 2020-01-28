#! /usr/local/anaconda3/bin/python
'''Generates artificial chap-seq reads'''

import argparse
import os
import sys
import numpy as np;
import random

from Bio import SeqIO

from afbio.generators import generator_doublesam
from afbio.numerictools import select_by_probability

parser = argparse.ArgumentParser(description='Generates artificial chap-seq reads');

parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to genome, fasta file");
parser.add_argument('--num_chimeras', nargs = '?', default=1000000, type = int, help = "Number of reads to generate");
parser.add_argument('--num_of_peaks', nargs = '?', default=300, type = int, help = "Number of peaks to generate");
parser.add_argument('--length', nargs = '?', default=300, type = int, help = "Length of the generated reads");
parser.add_argument('--min_strength', nargs = '?', default=3, type = int, help = "Minimum strength of the peaks");
parser.add_argument('--bs_width', nargs = '?', default=20, type = int, help = "Size of bindig sites");
args = parser.parse_args();



genome = SeqIO.to_dict(SeqIO.parse(args.genome, 'fasta'));
chrnames = list(genome.keys())
chrsizes = [len(genome[x]) for x in chrnames] 




def get_single(genome, chrnames, chrsizes, length):
    chrom = select_by_probability(chrnames, chrsizes);
    seqrecord = genome[chrom];
    ul = len(seqrecord) - length
    start = random.randint(0, ul);
    local = seqrecord[start:start+length]
    if(random.randint(0,1)):
        return str(local.seq.upper()), chrom, start, start+length, '+'
    else:
        return str(local.seq.reverse_complement().upper()), chrom, start, start+length, '-'

def get_doublet(genome, chrnames, chrsizes, length, mirnas, mlen = 13):
    mirid, mirseq = random.sample(mirnas.items(), 1)[0]
    mirseq = str(mirseq.seq.upper())
    l = random.randint(mlen, len(mirseq))
    mirseq = mirseq[:l]
    gseq, gchrom, gstart, gend, gstrand = get_single(genome, chrnames, chrsizes, length-l)
    seq = mirseq + gseq;
    
    return seq, gchrom, gstart, gend, gstrand, mirid
    

def convert(a):
    header = "|".join([ str(x) for x in a[1:]])
    return "@%s\n%s\n+\n%s" % (header, a[0], "F"*len(a[0]))
    
for _ in range(args.num_single):
    a = get_single(genome, chrnames, chrsizes, args.length)
    print(convert(a));
    
for _ in range(args.num_chimeras): 
    a = get_doublet(genome, chrnames, chrsizes, args.length, mirnas, mlen = 13)
    print(convert(a));
    
