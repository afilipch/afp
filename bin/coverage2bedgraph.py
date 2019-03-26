#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts coverage-per-position file into bedgrpah file'''

import argparse
import os
import sys


from pybedtools import BedTool
from collections import Counter

parser = argparse.ArgumentParser(description='Converts coverage-per-position file into bedgrpah file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the coverage file");
parser.add_argument('--trackopts', nargs = '?', default='', type = str, help = "Header for the (ucsc) track");
parser.add_argument('--chrom', nargs = '?', type = str, help = "If set, changes chromosome names to the provided one");
args = parser.parse_args();

if(args.trackopts):
    print(args.trackopts);


def flush(start, end, ident, chrom):
    if(ident and ident[1]):
        if(chrom):
            print("\t".join( [str(x) for x in (chrom, start-1, end-1, ident[1])] ))
        else:
            print("\t".join( [str(x) for x in (ident[0], start-1, end-1, ident[1])] ))

    
curstart = 0;
curident = None;
with open(args.path) as f:
    for l in f:
        a = l.strip().split("\t");
        start = int(a[1]);
        ident = (a[0], int(a[2]));
        if(ident != curident):
            flush(curstart, start, curident, args.chrom)
            curstart = start;
            curident = ident;
    else:
        flush(curstart, start, curident, args.chrom)
            
    
