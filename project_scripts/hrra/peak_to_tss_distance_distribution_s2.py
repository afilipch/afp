#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explore peak positions relative to the genes'''

import argparse
import os
import sys
#import copy
#from collections import defaultdict, Counter


import numpy as np;
from scipy.stats import rankdata
#import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Explore peak positions relative to the genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated peak regions");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Output destination for the plot");
parser.add_argument('--control', nargs = '?', required=True, type = str, help = "Path to the artificially generated peaks");
args = parser.parse_args();

def get_starts(path):
    peaks = BedTool(path);
    return [int(x.attrs['tss']) for x in peaks];

signal = get_starts(args.path)
control = get_starts(args.control)

brange = np.linspace(0, 1000, 11);
svals, cvals, xlabels = [], [], []
for s, e in zip(brange, brange[1:]):
    svals.append(len([x for x in signal if x>=s and x<e]))
    cvals.append(len([x for x in control if x>=s and x<e]))
    xlabels.append("<%d" % e);
    
svals.append(len([x for x in signal if x>=1000]))
cvals.append(len([x for x in control if x>=1000]))
norm = sum(svals)/sum(cvals)
cvals = [x*norm for x in cvals]
xlabels.append(">=1000");

print("distance\tnum of discovered peaks\tnum of random peaks");
for lsc in zip(xlabels, svals, cvals):
    print("%s\t%d\t%d" % lsc);


    
    
fig, ax = plt.subplots(figsize=(16, 9))
#plt.tight_layout(rect=[0.05, 0.05, 0.98, 0.98])
#ax.hist([signal, control], color = ['lightblue', 'coral'], normed = True, bins = 20, label = ['discovered peaks', 'random peaks'] );
ax.bar(brange, svals, 25, color='lightblue', label = 'discovered peaks')
ax.bar(brange+25, cvals, 25, color='coral', label = 'random peaks')
#plt.title(title)
plt.xticks(brange, xlabels)

#plt.xticks(brange, ["<%d" % x for x in lp[1:]] + [">%d" % lp[-1]], fontsize = 'xx-large')

plt.ylabel("number of peaks")
plt.xlabel("distance to the closest gene start")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(frameon=False, fontsize = 'xx-large')

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize('xx-large')

if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format=_format)
else:
    plt.show();
