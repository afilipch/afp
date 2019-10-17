#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores the evolution of total (averaged over all peaks) binding'''

import argparse
import os
import sys
from collections import defaultdict

import numpy as np;
import matplotlib.pyplot as plt
from pybedtools import BedTool, Interval



parser = argparse.ArgumentParser(description='Explores the evolution of total (averaged over all peaks) binding');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the consensus regions [regions/regions.annotated.gff]");
parser.add_argument('--zscore', nargs = '?', default=3, type = float, help = "Minimum required z-score");
parser.add_argument('--total', nargs = '?', const=True, default=False, type = str, help = "If set, the script reports total coverage, but not averaged");
parser.add_argument('--labels', nargs = '+', default=[], type = str, help = "Labels for the time series, for example: 0h 30m 2h");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");
args = parser.parse_args()

experiment2peaks = defaultdict(list);
for interval in BedTool(args.path):
    for num, (zscore, maxcov) in enumerate(zip(interval.attrs['zscores'].split(","), interval.attrs['maxcov'].split(","))):
        if(zscore != 'None' and maxcov != 'None'):
            #print(interval.attrs['zscores'])
            if(float(zscore) >= args.zscore):
                experiment2peaks[num].append(float(maxcov));

weights = [];
numpeaks = []


for key in sorted(experiment2peaks.keys()):
    numpeaks.append(len(experiment2peaks[key]));
    if(args.total):
        weights.append(sum(experiment2peaks[key]));
    else:
        weights.append(np.mean(experiment2peaks[key]));
    
for lw in zip(args.labels, weights):
    print("%s\t%1.2f" % lw)
                

### PLOTTING ###

fig, ax = plt.subplots(figsize=(16,9))
plt.xlabel('Time', fontsize='x-large')
if(args.total):
    plt.ylabel('Total peak height', fontsize='x-large')
else:
    plt.ylabel('Averaged peak height', fontsize='x-large')
    
plt.tick_params(axis='both', labelsize='x-large', top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

#set bins and boundaries
labels = args.labels
boundaries = range(0, len(labels));
bins = range(0, len(labels)+1);

#plot real and control binding pattern
plt.hist(boundaries, weights=weights, bins=bins, align='right', rwidth=0.6, color='skyblue')

#set xlabels
plt.xticks(range(1, len(labels)+1));
ax.set_xticklabels(labels, rotation=90)
ax.tick_params(axis='both', labelsize='x-large')


distance =  np.mean(weights)*0.05
for boundary, weight, text in zip(boundaries, weights, numpeaks):
    ax.text(boundary + 0.9, weight + distance, text, fontsize='xx-large');


if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format)
else:
    plt.show()
