#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Narrow down peak areas for better motif discovery'''

import argparse
import os
import sys
from collections import defaultdict, Counter


import numpy as np;
from pybedtools import BedTool, Interval
from afbio.generators import get_only_files



parser = argparse.ArgumentParser(description='Creates table with all the information regarding GLXR camp project');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with bed and/or gff files");
parser.add_argument('--mincov', nargs = '?', default=1, type = float, help = "Minimum allowed coverage")
parser.add_argument('--flank', nargs = '?', default=30, type = int, help = "Flank length of the peak area around its highest point")
parser.add_argument('--outdir', required=True, nargs = '?', type = str, help = "Path to the output directory")
args = parser.parse_args();

def check_interval(interval, mincov):
        return all( [float(x)>mincov for x in interval.attrs['topcoverage'].split(",")] );


for path in get_only_files(args.path):
    if(path.endswith('gff') or path.endswith('bed')):
        bedtool = BedTool(path)
        if(len(bedtool)):
            with open(os.path.join(args.outdir, os.path.basename(path)), 'w') as f:
                for interval in bedtool:
                    if(check_interval(interval, args.mincov)):
                        center = int(interval.name)
                        f.write(str(Interval(interval.chrom, center-args.flank, center+args.flank,  name = interval.name, strand = interval.strand, score=interval.attrs['topcoverage'])))
        
        
