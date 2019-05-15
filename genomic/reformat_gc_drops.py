#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Reformats gc drops into more convenient format'''

import argparse
import os
import sys
import pandas as pd;
from pybedtools import BedTool, Interval



parser = argparse.ArgumentParser(description='Explores relation between GC drops and binding peaks');
parser.add_argument('--gctrack', nargs = '?', required=True, type = str, help = "Path to the GC track");
parser.add_argument('--gcdrops', nargs = '?', required=True, type = str, help = "Path to the GC drops");
parser.add_argument('--maxdrop', nargs = '?', default=0.35, type = float, help = "Maximum allowed drop");
args = parser.parse_args();

gc_content = pd.read_csv(args.gctrack, sep="\t" , names = ["chr", "position", "gc"]).gc.values


#get gc drops
gcdrops = [];
with open(args.gcdrops) as f:
    current_drop = []
    for l in f:
        a = l.strip().split("\t");
        if(not current_drop and a[3] == 'max'):
            current_drop.append(int(a[1]));
        if(current_drop and a[3] == 'min'):
            current_drop.append(int(a[1]));
            gcdrops.append(tuple(current_drop));
            current_drop = [];
            
            
gcdrops_intervals = []             
for c, (start, end) in enumerate(gcdrops, start=1):
    score = min(gc_content[start: end])
    if(score <= args.maxdrop):
        sys.stdout.write(str(Interval('chr1', start, end, name = "drop_%d" % c, score = "%1.5f"  % score, strand = '+')));

