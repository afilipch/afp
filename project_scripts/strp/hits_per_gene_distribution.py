#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Draws a distribution of a number of bacterial hits per human membrane-gene'''

import argparse
import os
import sys
from collections import defaultdict
import copy


import numpy as np;
import matplotlib.pyplot as plt;




parser = argparse.ArgumentParser(description='Draws a distribution of a number of bacterial hits per human membrane-gene');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the processed blast hits, tsv format");
parser.add_argument('--query', nargs = '?', type = str, help = "Path to the all query sequences, fasta format");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Path to the output");
args = parser.parse_args()

gene2hits = defaultdict(int);
with open(args.path) as f:
    for l in f:
        query, hit = l.strip().split("\t")[:2]
        gene = "_".join(query.split(":")[:2])
        gene2hits[gene] += 1;
        
#print(len(gene2hits))


sys.exit(); 

motif_names = set([x[0] for x in BedTool(args.fimo)])
with_motif, without_motif = [], []




for peak in BedTool(args.path):
    intensity = np.mean([float(x) for x in peak.score.split(",")])
    if(peak.name in motif_names):
        with_motif.append(intensity);
    else:
        without_motif.append(intensity);
        
data = with_motif, without_motif
        
#print(len(motif_names))
#print(len(with_motif))

        
#sys.exit()


#############################################################################################################################
### DRAWING SECTION ###
fontsize = 24
linewidth = 3;
colors = ('lightblue', 'coral')
labels = ['with_motif\n(%s)' % len(data[0]), 'without_motif\n(%s)' % len(data[1])]
positions = [0.5,1]

fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

for datum, color, label, pos in zip(data, colors, labels, positions):
    boxplot = ax.boxplot(datum, positions = [pos], notch = True, patch_artist=True)
    for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(boxplot[item], color=color)
    plt.setp(boxplot["boxes"], facecolor=color)
    plt.setp(boxplot["fliers"], markeredgecolor=color)
    

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylabel('Averaged peak intensity', fontsize=fontsize)
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.set_xticklabels(labels)


if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format);
else:
    plt.show()
