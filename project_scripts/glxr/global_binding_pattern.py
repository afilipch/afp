#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Draws global binding patterns for multiple timepoints'''

import argparse
import os
import sys
#import copy
from collections import defaultdict, Counter


import numpy as np;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;
from afbio.generators import get_only_files
from afbio.sequencetools import coverage2dict

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Draws global binding patterns for various experiments');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the coverage folder");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Output destination for the plot");
parser.add_argument('--log', nargs = '?', default=False, const=True, type = bool, help = "If set, log2 transformation is apllied to the coverage");
parser.add_argument('--labels', nargs = '?', required=False, type = str, help = "Path to the file with peak labels");
parser.add_argument('--regions', nargs = '?', required=False, type = str, help = "Path to the annotated consensus regions");
parser.add_argument('--min_intensity', nargs = '?', default=3.0, type = float, help = "Minimum required peak intensity to be reported");
args = parser.parse_args();

### Labels
labels = [];

if(args.labels):
    with open(args.labels) as f:
        next(f)
        for l in f:
            a = l.strip().split("\t");
            name = a[6]
            if(name.startswith("NCgl")):
                name = a[5]
            #name.replace(" ", "\n")
            labels.append((name, int(a[2])))
elif(args.regions):
    for interval in BedTool(args.regions):
        if(interval.attrs['gtype'] == 'upstream'):
            score = max([float(x) for x in interval.attrs['topcoverage'].split(',')])
            if(score >= args.min_intensity):
                labels.append(( interval.attrs['genesymbol'], int(interval.name) ))
    #print(len(labels))






NAME_ORDER = ['glu_wt', 'glu_ko_cyab', 'ace_glu_wt', 'ace_glu_ko_cyab']


files = [x for x in get_only_files(args.path) if "normalized" in x]
name2files = defaultdict(list)
for f in files:
    name = "_".join(os.path.basename(f).split("_")[:-1])
    name2files[name].append(f)

size = len(name2files);
name2coverage = {}
for name, local_files in name2files.items():
    local_coverages = [list(coverage2dict(f).values())[0] for f in local_files]
    averaged_coverage = np.mean(local_coverages, axis=0)
    if(args.log):
        averaged_coverage = [np.log2(x+1) for x in averaged_coverage]
    name2coverage[name] = averaged_coverage
    length_coverage = len(averaged_coverage)
    


        



#############################################################################################################################
### DRAWING SECTION ###
fontsize = 28
colors = ('indigo', 'lawngreen', 'lightblue', 'coral')
dashline = 3;

positions = list(range(1, length_coverage+1));
fig, axes = plt.subplots(nrows=size, ncols = 1, sharex=False, sharey=True, figsize = (24, 7.5*size), frameon=False)
plt.tight_layout(rect=[0.05, 0.05, 0.96, 0.96], h_pad = 4)
#xticklocs, xticklabels = setticks(start, end)
ymax = max([max(x) for x in name2coverage.values()])*0.6
step = ymax*0.12
for ax, name, color in zip(axes, NAME_ORDER, colors):
    coverage = name2coverage[name]
    ax.plot(positions, coverage, color=color)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(positions[0], positions[-1]) 
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)
    ax.text(0.02, 0.95, name, verticalalignment='center', transform=ax.transAxes, fontsize=fontsize)
    ax.axhline(y=dashline, linestyle='--', color='dimgray', lw = 3)
    ax.set_xticklabels(['%1.1f' % (x/1000000) for x in ax.get_xticks()])
    
    
    adj_list = [x%5 for x in range(len(labels))]
    for (label, x), adj in zip(labels, adj_list):
        y = min(ymax, max(coverage)) + adj*step;
        y_peak = coverage[x]
        #print(y_peak)
        #ax.text(x, y, label, horizontalalignment='center', verticalalignment='center', fontsize=fontsize*0.75, style='italic')
        #ax.arrow(x,y,0,-2, width=2)
        if(y_peak > args.min_intensity):
            ax.annotate(label, xy=(x, y_peak), xytext=(x,y), arrowprops=dict(arrowstyle="->"), verticalalignment='center', fontsize=fontsize*0.75, style='italic')
        
    #ticks = [str(item[1]) for item in plt.xticks()]
    #print(ticks)
        
 
#plt.ylabel("enrichment factor");
fig.text(0.5, 0.01, 'genomic position (mB)', ha='center', fontsize=fontsize)

if(args.log):
    ytext = 'log2(normalized coverage)'
else:
    ytext = 'enrichment factor'
fig.text(0.005, 0.5, ytext, va='center', rotation='vertical', fontsize=fontsize)

if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format);
else:
    plt.show()
