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


parser = argparse.ArgumentParser(description='Draws a scatter plot peak intensities versus corresponding genes expression change');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the hrra all-in table");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Output destination for the plot");
##parser.add_argument('--control', nargs = '?', required=True, type = str, help = "Path to the artificially generated peaks");
args = parser.parse_args();

selected_names = ['hemH', 'hemE']

data = [];
with open(args.path) as f:
    m = next(f).strip().split("\t");
    for n, s in enumerate(m):
        print(n, s)
    for l in f:
        a = l.strip().split("\t")
        #if(a[9] != 'None' and a[21] != 'None' and float(a[6]) < 850):
            #data.append(( a[4], log(float(a[9]), 2), float(a[21]) ))
        if(a[8] != 'None' and a[9] != 'None' and float(a[6]) < 850):
            data.append(( a[4], log(float(a[8]), 2), log(float(a[9])/float(a[8]), 2) ))

name2coord = dict([ (x[0], x[1:]) for x in data ])


fontsize=24
linewidth = 5
xvals = np.array([x[1] for x in data])
yvals = np.array([x[2] for x in data])

fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.95])

ax.set_xlabel('Log2(peak intensity 0h)', fontsize=fontsize)
ax.set_ylabel('Log2(peak intensity 0.5h/0h)', fontsize=fontsize)    
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.scatter(xvals, yvals)
for name in selected_names:
    x, y = name2coord[name];
    ax.text(x,y,name)

if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format=_format)
else:
    plt.show();
