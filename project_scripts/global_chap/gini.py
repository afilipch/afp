#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''calculates gini distribution for peak scores'''

import argparse
import os
import sys
from os import listdir
from os.path import isfile
from collections import defaultdict

import numpy as np;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

from afbio.pybedtools_af import read_comments




parser = argparse.ArgumentParser(description='Explores shared targets among different dna-binding proteins');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the merged binding regions, gff files");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
#parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to NCBI annotation table");
#parser.add_argument('--maxsd', nargs = '?', default = 350, type = int, help = "Maximal allowed distance to the start gene");
#parser.add_argument('--strict', nargs = '?', default = False, const=True, type = bool, help = "If set, strict requirments for a gene to be counted as controlled by an RBP are applied");
args = parser.parse_args();



position2name = dict([ (x[0], x[1].split("_")[0]) for x in enumerate(read_comments(args.path)[0][2:].split(",")) ])




protein2intensities = defaultdict(list)
for interval in BedTool(args.path):
    lints = [float(x) for x in interval.attrs['topcoverage'].split(",")]
    locald = defaultdict(list)
    for p, val in enumerate(lints):
        if(val):
            locald[position2name[p]].append(val);
    for name, vals in locald.items():
        protein2intensities[name].append(np.mean(vals))
        
labels = []
values = []
for k, v in protein2intensities.items():
    gini = np.percentile(v, 75)/np.percentile(v, 25)
    labels.append(k)
    values.append(gini)



### DRAWING

fontsize = 24
linewidth = 3

x = np.arange(len(values))
barcolors = ['darkblue', 'lightblue', 'darkgray', 'lightgray']
barlabels = ['phage', 'non-phage', 'control phage', 'control non-phage']


fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.9])

ax.set_xlabel('Protein', fontsize=fontsize)
ax.set_ylabel('Gini_25 index', fontsize=fontsize)
ax.set_xticks(x)
ax.set_xticklabels(labels)    
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for axis in ['bottom','left','right']:
    ax.spines[axis].set_linewidth(linewidth)

ax.bar(x, values, 0.8)
plt.savefig(os.path.join(args.outdir, "gini.%s"  % args.format) , format = args.format)

    





    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
