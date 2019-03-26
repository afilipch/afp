#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts gene ids to gene symbol in a given file'''

import argparse
import os
import sys
from collections import defaultdict


from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Finds and explores differentially expressed genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the file with gene ids");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the gene annotation file");
parser.add_argument('--ncol', nargs = '?', required=True, type = int, help = "Position of the column with gene ids");
parser.add_argument('--header', nargs = '?', default=0, type = int, help = "Number of header lines");
args = parser.parse_args();

geneid2genesymbol = {};
for interval in BedTool(args.annotation):
    geneid2genesymbol[interval.attrs['geneid']] = interval.attrs['Name']
    
with open(args.path) as f:
    for _ in range(args.header):
        print(next(f).strip());
    for l in f:
        a = l.strip().split("\t");
        a[args.ncol] = geneid2genesymbol.get(a[args.ncol], a[args.ncol])
        print("\t".join(a));
#for interval in BedTool(args.path):
    #interval
    