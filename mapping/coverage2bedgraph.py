#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts coverage-per-position file into bedgrpah file'''

import argparse
import os
import sys


from pybedtools import BedTool
from collections import Counter
import numpy as np

from afbio.sequencetools import coverage2dict;

parser = argparse.ArgumentParser(description='Converts coverage-per-position file into bedgrpah file');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the coverage file");
parser.add_argument('--trackopts', nargs = '?', default='', type = str, help = "Header for the (ucsc) track");
parser.add_argument('--multiplier', nargs = '?', default=1.0, type = float, help = "Coverage multiplier, should be used only for the normalized tracks");
parser.add_argument('--convert', nargs = '?', default=False, const=True, type = bool, help = "If set, converts chromosome names to human ones");
parser.add_argument('--normalize', nargs = '?', default=False, const=True, type = bool, help = "If set, coverage is normalized to the mean value");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Output destination for the statistics plots");
args = parser.parse_args();


def flush(start, end, val, chr_count):
    print("\t".join( [str(x) for x in ('chr%d' % chr_count, start, end, val)] ))


def convert_coverage(coverage, chr_count):
    #print(len(coverage))
    curstart = 0;
    curval = coverage[0];

    for pos, val in enumerate(coverage[1:], start=1):
        if(val != curval):
            #print(val, curval, curstart)
            flush(curstart, pos, curval, chr_count)
            curstart = pos;
            curval = val;
    else:
        flush(curstart, pos, val, chr_count)
            
            



if(args.trackopts):
    print(args.trackopts);
    
    

coverage_dict_list = [coverage2dict(x) for x in args.path];
coverage_dict = {};
for chrom in coverage_dict_list[0].keys():
    temp = sum([x[chrom] for x in coverage_dict_list])
    if(args.normalize):
        temp = temp/np.mean(temp)
    coverage_dict[chrom] = [int(x) for x in temp*args.multiplier]
    
coverage_dict = [x[1] for x in sorted(coverage_dict.items(), key = lambda x: x[0])]
for chr_count, coverage in enumerate(coverage_dict, start=1):
    convert_coverage(coverage, chr_count)
    

#sys.exit();







#curstart = 0;
#curident = None;
#cur_chrom = None
#chrom_ident = None
#chrom_end = None

#with open(args.path) as f:
    #for l in f:
        #a = l.strip().split("\t");
        #start = int(a[1]);
        #ident = (a[0], round(float(a[2])*args.multiplier));
        #if(cur_chrom and  a[0] != cur_chrom):
            #flush(curstart, chrom_end, chrom_ident, args.convert, chr_count)
            #chr_count -= 1;
            #curstart = start;
            #curident = ident

        #elif(ident != curident):
            #flush(curstart, start, curident, args.convert, chr_count)
            #curstart = start;
            #curident = ident;
        #chrom_ident = ident;
        #chrom_end = start+1
        #cur_chrom = a[0]
    #else:
        #flush(curstart, start, curident, args.convert, chr_count)
            
    
