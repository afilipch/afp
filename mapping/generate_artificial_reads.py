#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Generates artificial paired end reads from a given genome'''

import argparse
import os
import sys
import numpy as np;
import random

from Bio import SeqIO;
from Bio.Seq import reverse_complement;

from afbio.generators import generator_doublesam
from afbio.numerictools import select_by_probability

parser = argparse.ArgumentParser(description='Generates artificial reads from a given genome');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome, fasta file");
parser.add_argument('--output', nargs = '?', required=True, type = str, help = "Base name for the output files. If set '--output test', then 'test.r1.fastq and test.r2.fastq are created'");
parser.add_argument('-fl', '--fragment_length', nargs = '?', default=200, type = int, help = "Length of the generated fragments");
parser.add_argument('--length', nargs = '?', default=60, type = int, help = "Length of the generated paired-end reads");
parser.add_argument('--num', nargs = '?', default=100000, type = int, help = "Number of reads to generate");
args = parser.parse_args();

genome = SeqIO.to_dict(SeqIO.parse(args.path, 'fasta'));

chrnames = list(genome.keys())
chrsizes = [len(genome[x]) for x in chrnames] 

def get_single(genome, chrnames, chrsizes, length):
    chrom = select_by_probability(chrnames, chrsizes);
    seqrecord = genome[chrom];
    ul = len(seqrecord) - length-1
    start = random.randint(1, ul);
    end = start+length
    local = str(seqrecord[start:end].seq.upper())
    lflank = str(seqrecord[start-1:start+1].seq.upper())
    rflank = str(seqrecord[end-1: end+1].seq.upper())
    if(random.randint(0,1)):
        return local, chrom, start, end, '+', lflank, rflank
    else:
        return reverse_complement(local), chrom, start, end, '-', reverse_complement(rflank), reverse_complement(lflank)
    
    
    
    
def convert(a, count, length):
    header = "|".join([ str(x) for x in a[1:]])
    return "@%d|%s\n%s\n+\n%s\n" % (count, header, reverse_complement(a[0][-length:]), "F"*length), "@%d|%s\n%s\n+\n%s\n" % (count, header, a[0][:length], "F"*length)

    
    
    
    
with open(args.output+'.r1.fastq', 'w') as r1, open(args.output+'.r2.fastq', 'w') as r2:
    for count in range(args.num):
        a = get_single(genome, chrnames, chrsizes, args.fragment_length);
        s1, s2 = convert(a, count+1, args.length);
        r1.write(s1);
        r2.write(s2);
        
    
    