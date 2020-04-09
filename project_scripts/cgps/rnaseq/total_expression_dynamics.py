#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Draws the chart of expression dynamics over time for phage VS non-phage genes'''

import argparse
import os
import sys
from collections import defaultdict
import copy


import numpy as np;
from pybedtools import BedTool, Interval
from itertools import combinations, permutations
import matplotlib.pyplot as plt;




parser = argparse.ArgumentParser(description='Draws the chart of expression dynamics over time for phage VS non-phage genes');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the expression transcripts (first non-phage, second phage), gff format");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path to the output");
args = parser.parse_args()

data = [];
for path in args.path:
    temp = [];
    with open(path) as f:
        xlabels = next(f).strip().split("=")[1].split(",")
    temp = [0 for _ in xlabels];
    for interval in BedTool(path):
        for c, el in enumerate(interval.attrs['expression'].split(':')):
            temp[c] += np.mean([float(x) for x in el.split(",")])
    data.append(temp);
        #data.append(  np.array([float(x) for x in next(f).strip().split("=")[1].split(",")])/1000000  )
        





#############################################################################################################################
### DRAWING SECTION ###
fontsize = 24
linewidth = 3;
colors = ('lightblue', 'coral')
labels = ['non-phage', 'phage']


fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

for datum, color, label in zip(data, colors, labels):
    ax.plot(datum, label=label, color=color, linewidth=linewidth)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("Time", fontsize=fontsize)
ax.set_ylabel('Total expression [fraction]', fontsize=fontsize)
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.set_xticklabels([''] + xlabels)

fig.legend(fontsize=fontsize, frameon=False, loc = (0.7, 0.8))

if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format);
else:
    plt.show()
