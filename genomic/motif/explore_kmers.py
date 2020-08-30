#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores kmers frequences content in the provided sequences'''

import argparse
import os
import sys
from collections import defaultdict

import numpy as np;
from Bio import SeqIO

from afbio.sequencetools import sliding_window



parser = argparse.ArgumentParser(description='Explores kmers frequences content in the provided sequences');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the sequences, fasta format");
parser.add_argument('--kmer', nargs = 2, default=[8,8], type = int, help = "Min and max length of the kmer, both are included");
parser.add_argument('--top', nargs = '?', default=100, type = int, help = "Output N the most frequent kmers");
parser.add_argument('--reverse', nargs = '?', default=False, const=True, type = bool, help = "If set, kmers on reverse strand are also counted");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();

kmer_range = range(args.kmer[0], args.kmer[1] + 1)

for size in kmer_range:
    kmers_count = defaultdict(int);
    norma = 0;
    with open(os.path.join(args.outdir, 'kmer_%d.tsv' % size), 'w') as f:
        for seqrecord in SeqIO.parse(args.path, 'fasta'):
            seqlength = len(seqrecord)
            norma += 1;
            curset = set()
            for kmer in sliding_window(seqrecord.seq, size):
                curset.add(kmer)
            if(args.reverse):
                for kmer in sliding_window(seqrecord.seq.reverse_complement(), size):
                    curset.add(kmer)
            for el in curset:
                kmers_count[el] += 1;
        meanfr = (seqlength - size + 1)/(4**size)
        if(args.reverse):
            meanfr *= 2;
                
        kmers_count = [(x[0], x[1]/norma) for x in kmers_count.items()]
        kmers_count.sort(key = lambda x: x[1], reverse = True)
        for kmer, fraction in kmers_count[:args.top]:
            f.write("%s\t%f\t%1.2f\n" % (''.join(kmer), fraction, fraction/meanfr))
    #print(meanfr)


            
    


