#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explore peak positions relative to the genes'''

import argparse
import os
import sys
#import copy
#from collections import defaultdict, Counter


import numpy as np;
from math import log;
#from scipy.stats import rankdata
#import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Draws a scatter plot peak TSS distances versus corresponding genes expression change');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the hrra all-in table");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Output destination for the plot");
args = parser.parse_args();



data = [];
with open(args.path) as f:
    m = next(f).strip().split("\t");
    for n, s in enumerate(m):
        print(n, s)
    for l in f:
        a = l.strip().split("\t")
        expr = float(a[21])
        tss = float(a[6])
        if(abs(expr) < 10):
            data.append((tss, expr))




fontsize=24
linewidth = 5
xvals = np.array([x[0] for x in data])
yvals = np.array([x[1] for x in data])

fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.95])

ax.set_xlabel('TSS distance', fontsize=fontsize)
ax.set_ylabel('Log2(peak intensity 0.5h KO/WT)', fontsize=fontsize)    
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.scatter(xvals, yvals)


if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format=_format)
else:
    plt.show();
