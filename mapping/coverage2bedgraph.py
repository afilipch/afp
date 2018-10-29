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
parser.add_argument('--multiplier', nargs = '?', default=1.0, type = float, help = "Coverage multiplier, should be used only for the normalized tracks");
parser.add_argument('--convert', nargs = '?', default=False, const=True, type = bool, help = "If set, converts chromosome names to human ones");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Output destination for the statistics plots");
args = parser.parse_args();

if(args.trackopts):
    print(args.trackopts);


def flush(start, end, ident, convert):
    if(ident and ident[1]):
        if(convert):
            print("\t".join( [str(x) for x in ('chr1', start, end, ident[1])] ))
        else:
            print("\t".join( [str(x) for x in (ident[0], start, end, ident[1])] ))

    
curstart = 0;
curident = None;
with open(args.path) as f:
    for l in f:
        a = l.strip().split("\t");
        start = int(a[1]);
        ident = (a[0], round(float(a[2])*args.multiplier));
        if(ident != curident):
            flush(curstart, start, curident, args.convert)
            curstart = start;
            curident = ident;
    else:
        flush(curstart, start, curident, args.convert)
            
    
