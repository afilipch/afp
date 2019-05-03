#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Calculates smoothed (averaged) GC content along the provided genome'''

import argparse
import sys
import os
import numpy as np;
from Bio import SeqIO
from collections import Counter;

from afbio.sequencetools import sliding_window

parser = argparse.ArgumentParser(description='Calculates smoothed (averaged) GC content along the provided genome');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome, fasta format");
parser.add_argument('--flen', nargs = '?', default=25, type = int, help = "Flank length to calculate an average GC content (window length is equal to flen*2+1)");
args = parser.parse_args();


wlen = args.flen*2+1
args = parser.parse_args();

def get_flanks(window):
    for c in range(int(len(window)/2)):
        yield window[:c*2+1];
        
def get_flanks_backward(window):
    for c in range(int(len(window)/2)):
        yield window[c*2+2:];        

def get_content(window):
    counter = Counter(window);
    return (counter['G'] + counter['C'])/float(len(window));

for seqrec in SeqIO.parse(args.path, 'fasta'):
    chrom = seqrec.id
    seq = seqrec.seq.upper();
    position = 0;
    for count, window in enumerate(sliding_window(seq, wlen)):
        if(count == 0):
            for flank in get_flanks(window):
                sys.stdout.write("%s\t%d\t%1.5f\n" % (chrom, position, get_content(flank)))
                position += 1;
        sys.stdout.write("%s\t%d\t%1.5f\n" % (chrom, position, get_content(window)))
        position += 1;
    else:
        #sys.stderr.write("%s\n" % (window,))
        for flank in get_flanks_backward(window):
            sys.stdout.write("%s\t%d\t%1.5f\n" % (chrom, position, get_content(flank)))
            position += 1;
        
