#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Reformats AT extrema into more convenient format'''

import argparse
import os
import sys
import pandas as pd;
from pybedtools import BedTool, Interval

from collections import defaultdict
import numpy as np

from afbio.sequencetools import coverage2dict;



parser = argparse.ArgumentParser(description='Reformats AT extrema into more convenient format');
parser.add_argument('--attrack', nargs = '?', required=True, type = str, help = "Path to the AT track");
parser.add_argument('--atdrops', nargs = '?', required=True, type = str, help = "Path to the AT peaks");
parser.add_argument('--minat', nargs = '?', default=0.5, type = float, help = "Manimum allowed maximum AT content inside a peak");
args = parser.parse_args();

at_content_dict = coverage2dict(args.attrack)


#get gc drops
gcdrops = defaultdict(list);
with open(args.atdrops) as f:
    current_drop = []
    for l in f:
        a = l.strip().split("\t");
        if(not current_drop and a[3] == 'max'):
            current_drop.append(int(a[1]));
        if(current_drop and a[3] == 'min'):
            current_drop.append(int(a[1]));
            gcdrops[a[0]].append(tuple(current_drop));
            current_drop = [];
            
                         
for chrom, positions in gcdrops.items():
    for c, (start, end) in enumerate(positions, start=1):
        score = max(at_content_dict[chrom][start: end])
        if(score >= args.minat):
            sys.stdout.write(str(Interval(chrom, start, end, name = "drop_%s_%d" % (chrom, c), score = "%1.3f"  % score, strand = '+')));

