#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explore the correlation of samples (along with replicates) for the gene expression'''

import argparse
import os
import sys
from itertools import combinations;
#import copy
#from collections import defaultdict, Counter


import numpy as np;
from scipy.stats import pearsonr,spearmanr
#import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Explore evolution of the peaks over time');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the differential table, tsv file");
parser.add_argument('--mincov', nargs = '?', default=10, type = float, help = "Minimum required coverage [TPM]");
parser.add_argument('--log', nargs = '?', default=False, const = True, type = bool, help = "If set, log transformation is applied");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Path to the plot");
args = parser.parse_args();


if(args.log):
    log_string = ' log(n+1)'
else:
    log_string = ''

### Read input data
START = 4
with open(args.path) as f:
    labels = [x for x in next(f).strip().split("\t")[START:] if "|" not in x]
    stop = START + len(labels)
    expr_list = [[] for i in range(len(labels))];
    #print(expr_list)
    for l in f:
        a = [x for x in l.strip().split("\t")[START:stop]]
        #print(a)
        for c, el in enumerate(a):
            expr_list[c].append([ float(x) for x in el.split(",")])
            

print("sample\tbatch\tbatch\tPearson R%s" % log_string)
expr_list = [np.array(x) for x in expr_list]
for sample, label in zip(expr_list, labels):
    for c1, c2 in combinations(range(sample.shape[1]), 2):
        a1, a2 = [], []
        for el1, el2 in zip(sample[:,c1], sample[:,c2]):
            if(el1>=args.mincov or el2>=args.mincov):
                a1.append(el1)
                a2.append(el2)
                
        if(args.log):
            a1 = [np.log(x+1) for x in a1]
            a2 = [np.log(x+1) for x in a2]
            
        pc, pval = pearsonr(a1, a2)
        #pc, pval = spearmanr(a1, a2)
        #print("%s\t%d\t%d\t%1.3f" % (label, c1+1, c2+1, pc))
    #print()
print()







######################################################################################
### Find correlation coefficients

inter_samples = (np.array([np.mean(x, axis=1) for x in expr_list])).transpose()
print("sample1\tsample2\tPearson R%s" % log_string)
cmatrix = np.zeros( (inter_samples.shape[1], inter_samples.shape[1]));
for i in range(cmatrix.shape[0]):
    cmatrix[i] = 1
    


for i, j in combinations(range(inter_samples.shape[1]), 2):
    with_signal = np.array( [(x[0], x[1]) for x in zip(inter_samples[:,i], inter_samples[:,j]) if max(x)>args.mincov])
    if(args.log):
        with_signal = np.array([(np.log(x[0]+1), np.log(x[1]+1)) for x in with_signal])
        
    pc, pval = pearsonr(with_signal[:,0], with_signal[:,1])
    cmatrix[i,j] = pc
    cmatrix[j,i] = pc
    print("%s\t%s\t%1.3f" % (labels[i], labels[j], pc) )



######################################################################################
### Draw a heatmap

timestamps = labels
cmap="RdPu"

fig, ax = plt.subplots()
plt.tight_layout(rect=[0.05, 0, 1, 0.90])
im = ax.imshow(cmatrix, cmap=cmap, vmax=1)
cbar = ax.figure.colorbar(im, ax=ax, cmap=cmap)
cbar.ax.set_ylabel("Pearson correlation%s" % log_string, rotation=-90, va="bottom")

ax.set_xticks(np.arange(len(timestamps)))
ax.set_yticks(np.arange(len(timestamps)))
ax.set_xticklabels(timestamps, rotation=-90)
ax.set_yticklabels(timestamps)
ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

for edge, spine in ax.spines.items():
    spine.set_visible(False)

ax.set_xticks(np.arange(cmatrix.shape[1]+1)-.5, minor=True)
ax.set_yticks(np.arange(cmatrix.shape[0]+1)-.5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
ax.tick_params(which="minor", bottom=False, left=False)


if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format)
else:
    plt.show()
