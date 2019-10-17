#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Draws global binding patterns for multiple timepoints'''

import argparse
import os
import sys
#import copy
#from collections import defaultdict, Counter


import numpy as np;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;
import pandas as pd
from math import log

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Draws global binding patterns for multiple timepoints');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the normalized coverage, bed format");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Output destination for the plot");
parser.add_argument('--log', nargs = '?', default=False, const=True, type = bool, help = "If set, log2 transformation is apllied to the coverage");
args = parser.parse_args();

size = len(args.path);

if(not args.log):
    coverage_set = [pd.read_csv(path, sep="\t" , names = ["chr", "postion", "coverage"]).coverage.values for path in args.path]
else:
    coverage_set = []
    for path in args.path:
        coverage = [log(x+1,2) for x in pd.read_csv(path, sep="\t" , names = ["chr", "postion", "coverage"]).coverage.values]
        coverage_set.append(np.array(coverage));



#print(max(coverage_set[1]));
#sys.exit()

positions = list(range(1, len(coverage_set[0])+1));


fig, axes = plt.subplots(nrows=size, ncols = 1, sharex=False, sharey=True, figsize = (16, 5*size), frameon=False)
plt.tight_layout(rect=[0.05, 0.05, 0.98, 0.98], h_pad = 4)
#xticklocs, xticklabels = setticks(start, end)
for ax, coverage in zip(axes, coverage_set):
    ax.plot(positions, coverage, 'b')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(positions[0], positions[-1])
    #ax.spines['bottom'].set_visible(False)
    #ax.set_xticks([]) 
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()   
    ax.tick_params(axis='both', labelsize='xx-large')
    
fig.text(0.5, 0.01, 'genomic position (nt)', ha='center', fontsize='xx-large')

if(args.log):
    ytext = 'log2(normalized coverage)'
else:
    ytext = 'normalized coverage'
fig.text(0.005, 0.5, ytext, va='center', rotation='vertical', fontsize='xx-large')

if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format);
else:
    plt.show()
