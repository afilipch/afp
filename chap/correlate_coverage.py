#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Correlates chap/chip genomic coverage among the given samples'''

import argparse
import os
import sys
from scipy.stats import pearsonr
from itertools import combinations
import numpy as np;
import pandas as pd;
import matplotlib.pyplot as plt;


parser = argparse.ArgumentParser(description='Correlates chap/chip genomic coverage among the given samples');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the coverage tracks");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Path for the output plot");
args = parser.parse_args();

### Read the input
coverage_list = [pd.read_csv(x, sep="\t" , names = ["chr", "postion", "coverage"]).coverage.values for x in args.path]
names = [ ".".join(os.path.basename(x).split('.')[:-1]) for x in args.path]



######################################################################################
### Find correlation coefficients

cmatrix = np.zeros((len(coverage_list), len(coverage_list)));
for i in range(cmatrix.shape[0]):
    cmatrix[i] = 1

for i, j in combinations(range(cmatrix.shape[0]), 2):
    pc = pearsonr(coverage_list[i], coverage_list[j])[0]
    cmatrix[i,j] = pc
    cmatrix[j,i] = pc
    sys.stderr.write("correlation between pattern %d and %d is equal %1.3f\n" % (i+1, j+1, pc))
    
for r in cmatrix:
    print("\t".join(['%1.3f' % x for x in r]))


######################################################################################
### Draw a heatmap

timestamps = names
cmap="GnBu"

fig, ax = plt.subplots()
im = ax.imshow(cmatrix, cmap=cmap, vmax=1)
cbar = ax.figure.colorbar(im, ax=ax, cmap=cmap)
cbar.ax.set_ylabel('Pearson correlation', rotation=-90, va="bottom")

ax.set_xticks(np.arange(len(timestamps)))
ax.set_yticks(np.arange(len(timestamps)))
ax.set_xticklabels(timestamps)
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