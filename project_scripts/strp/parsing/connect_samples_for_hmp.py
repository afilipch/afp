#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Connects samples names for identical proteins'''

import argparse
import os
import sys

from Bio import SeqIO


parser = argparse.ArgumentParser(description='Connects samples names for identical proteins');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to protein fasta file");
args = parser.parse_args()



        
#print(len(selected))
for seqrecord in SeqIO.parse(args.path, 'fasta'):
    l = seqrecord.description.split(" ")
    if(len(l) == 1):
        print(l[0])
    else:
        print("\t".join(l[:1] + l[1].split(",")))

