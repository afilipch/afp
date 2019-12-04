#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explore evolution of the peaks over time'''

import argparse
import os
import sys
from itertools import combinations;
#import copy
#from collections import defaultdict, Counter


import numpy as np;
from scipy.stats import pearsonr
#import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Explore evolution of the peaks over time');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the consensus regions, gff file");
parser.add_argument('--names', nargs = '+', required=True, type = str, help = "Names of the provided samples, must be of the same order as values in \'zscores\' attribute of the provided consensus regions");
parser.add_argument('--min-zscore', nargs = '?', default=1, type = float, help = "Minimum required zscore for both peaks to calculate correlation between them");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Path to the output directory for the plots");
parser.add_argument('--selection', nargs = '+', default=None, type = int, help = "If set only the selected (by number starting from one) sets of peak will be correlated");
#parser.add_argument('--genes', nargs = '?', default='', type = str, help = "Path to the selected genes to check for correlation, tsv format");
args = parser.parse_args();

### Read input data

LOWEST_ZSCORE = -10;
COV_THRESHOLD = 0;
min_zscore = max(args.min_zscore, LOWEST_ZSCORE); 
regions = BedTool(args.path);


def get_local_coverages(region, min_zscore):
    temp_covs = [float(x) if x != 'None' else -1 for x in region.attrs['maxcov'].split(",")]
    temp_zscores = [float(x) if x != 'None' else LOWEST_ZSCORE for x in region.attrs['zscores'].split(",")]
    covs = [x[0] if x[1] >= min_zscore else -1 for x in zip(temp_covs, temp_zscores)]
    #print(covs)
    #print(region)
    #print("________________________________________________________________________________________________________")
    return covs;


maxcovs = np.array([get_local_coverages(region, min_zscore) for region in regions])
if(args.selection):
    selection = [x-1 for x in args.selection]
    maxcovs = maxcovs[:,selection]
    

#print(maxcovs[:,0])






######################################################################################
### Find correlation coefficients

print("\n")
cmatrix = np.zeros( (maxcovs.shape[1], maxcovs.shape[1]));
for i in range(cmatrix.shape[0]):
    cmatrix[i] = 1

for i, j in combinations(range(maxcovs.shape[1]), 2):
    with_signal = np.array([(x[0], x[1]) for x in zip(maxcovs[:,i], maxcovs[:,j]) if min(x)>COV_THRESHOLD])
    pc = pearsonr(with_signal[:,0], with_signal[:,1])[0]
    cmatrix[i,j] = pc
    cmatrix[j,i] = pc
    print("correlation between pattern %d and %d is equal %.3f. Num of compared peaks: %d" % (i+1, j+1, pc, len(with_signal)))



######################################################################################
### Draw a heatmap

timestamps = args.names
cmap="RdPu"

fig, ax = plt.subplots()
plt.tight_layout(rect=[0.2, 0, 1, 0.7])
im = ax.imshow(cmatrix, cmap=cmap, vmax=1)
cbar = ax.figure.colorbar(im, ax=ax, cmap=cmap)
cbar.ax.set_ylabel('Pearson correlation', rotation=-90, va="bottom")

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


    
