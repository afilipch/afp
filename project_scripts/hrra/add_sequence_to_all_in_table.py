#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Adds genomic sequences to the hrra all in table'''

import argparse
import os
import sys

import numpy as np
from pybedtools import BedTool, Interval
from collections import defaultdict, Counter
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Adds genomic sequences to the hrra all in table');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the hrra all-in table");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome");
args = parser.parse_args()

genome = str(next(SeqIO.parse(args.genome, 'fasta')).seq)



INS_PLACE = 3
taken = set()
with open(args.path) as f:
    header = next(f).strip().split("\t")
    header.insert(INS_PLACE, 'Sequence')
    print("\t".join(header))
    for l in f:
        a = l.strip().split("\t")
        start, stop = int(a[1]), int(a[2])
        if((start, stop) not in taken):
            a.insert(INS_PLACE, genome[start:stop])
            taken.add((start, stop))
            print("\t".join(a))
        
    

            
            

