#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores relation between AT stretches number and binding peaks intensities'''

import argparse
import os
import sys
import numpy as np;
from scipy.stats import pearsonr
import pandas as pd;
from pybedtools import BedTool, Interval
from Bio import SeqIO




parser = argparse.ArgumentParser(description='Explores relation between AT stretches and binding peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks");
parser.add_argument('--atstretches', nargs = '?', required=True, type = str, help = "Path to the AT stretches, bed format");
#parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");
parser.add_argument('--minlength', nargs = '?', default=7, type = int, help = "Minimal allowed length of an AT stretch");

args = parser.parse_args();

regions = BedTool(args.path);
atstretches = BedTool( [ x for x in BedTool(args.atstretches) if int(x.score)<1 and len(x) >= args.minlength] )

pairs = []
for interval in regions.intersect(atstretches, c = True):
    pairs.append((float(interval.score), float(interval.attrs['topcoverage']), float(interval[-1])))
zscores = [x[0] for x in pairs];
topcoverages = [x[1] for x in pairs];
nums = [x[2] for x in pairs];


zlimits = [2, 5, 10, 25, 100, 1000]
zlimit2means = [];
climits = [1, 2.5, 5, 10, 50, 1000]
climit2means = [];

zcounts = [] 
ccounts = []

for l1, l2 in zip(zlimits, zlimits[1:]):
    zlimit2means.append(np.mean([x[2] for x in pairs if x[0] >= l1 and x[0] < l2]))
    zcounts.append(len([x for x in pairs if x[0] >= l1 and x[0] < l2]))
for l1, l2 in zip(climits, climits[1:]):
    climit2means.append(np.mean([x[2] for x in pairs if x[1] >= l1 and x[1] < l2]))
    ccounts.append(len([x for x in pairs if x[1] >= l1 and x[1] < l2]))
    
    
print(pearsonr(zscores, nums));
print(zlimit2means);
print(zcounts);
print(climit2means);
print(ccounts);

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(16,9))

barwidth = 0.9;
labels = ["%1.1f-%1.1f" % x for x in zip(climits, climits[1:])]
uplabels = ["n=%d" % x for x in ccounts]
bars = climit2means
positions = np.arange(len(bars));
plt.bar(positions, bars, width = barwidth, color = (0.3,0.1,0.4,0.6))
for pos, bar, label in zip(positions, bars, uplabels):
    plt.text(x = pos-barwidth/3 , y = bar+0.1, s = label, size = 18)

plt.ylabel('Number of AT stretches per peak')
plt.xlabel("Peak top coverage [averaged coverage units]")
plt.xticks(positions, labels, rotation=90)
plt.subplots_adjust(bottom= 0.2, top = 0.98)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False) 
 
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize('x-large') 
 
if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format)
else:
    plt.show()
