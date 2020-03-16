#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Selects human protein sequences based on the provided tsv file'''

import argparse
import os
import sys

from Bio import SeqIO


parser = argparse.ArgumentParser(description='Selects human protein sequences based on the provided tsv file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to protein fasta file");
parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Path to the table with proteins from human genome atlas");
args = parser.parse_args()


selected = set();
with open(args.table) as f:
    header = next(f).strip().split("\t");
    uniprot = header.index('Uniprot')
    for l in f:
        selected.add(l.strip().split("\t")[uniprot])
        
#print(len(selected))
for seqrecord in SeqIO.parse(args.path, 'fasta'):
    name = seqrecord.name.split("|")[1]
    if(name in selected):
        print(">%s\n%s" % (seqrecord.name, seqrecord.seq) )
