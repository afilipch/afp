#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Evaluates gene expression values based on coverage'''

import argparse
import os
import sys
from collections import defaultdict, namedtuple
import math
from itertools import combinations

import numpy as np;
import pandas as pd;
import scipy 
import matplotlib.pyplot as plt;
from scipy.stats.stats import pearsonr, spearmanr;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Evaluates gene expression values based on coverage');
#parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the sam file");
parser.add_argument('--first', nargs = '+', required=True, type = str, help = "Path to the genes' expression in the first condition, tsv format");
parser.add_argument('--second', nargs = '+', required=True, type = str, help = "Path to the genes' expression in the second condition, tsv format");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Path to the output directory for the plots");
parser.add_argument('--corrtype', nargs = '?', default='pearson', choices = ['spearman', 'pearson'], type = str, help = "Type of correlation, spearman or pearson");

args = parser.parse_args();

#gene2counts = defaultdict(lambda: defaultdict(list));

def get_tpms(multipath):
    vecs = [];
    for path in multipath:
        vec = [];
        genes = [];
        with open(path) as f:
            for l in f:
                a = l.strip().split("\t");
                vec.append(float(a[1]))
                genes.append(a[0])
        vecs.append(np.array(vec));
    vecs = np.array(vecs);
    vecs = vecs.transpose()
    return genes, vecs
        

genes, m_first = get_tpms(args.first);
genes, m_second = get_tpms(args.second);
m_total = np.concatenate((m_first, m_second), axis = 1)


###Estimate intra/inter -sample correlations for gene expression 
if(args.corrtype == 'spearman'):
    corfun = spearmanr
if(args.corrtype == 'pearson'):
    corfun = pearsonr



def get_correlations(array):
    correlation = np.zeros((array.shape[1], array.shape[1]))
    for i in range(correlation.shape[0]):
        correlation[i, i] = 1
        
    for i1, i2 in combinations(range(array.shape[1]), 2):
        pc_coeff = corfun(array[:,i1], array[:,i2])[0]
        correlation[i1, i2] = pc_coeff;
        correlation[i2, i1] = pc_coeff;
        #print(array.shape)
        
    return correlation;
        
correlation = get_correlations(m_total);
#print(correlation)
#sys.exit()





######################################################################################
### Draw a heatmap
def get_stamps(multipath):
    return ["%s\n%s" % tuple(x.split(".")[:2]) for x in multipath];

timestamps = get_stamps(args.first) + get_stamps(args.second)
fig, ax = plt.subplots()
im = ax.imshow(correlation, cmap="RdPu", vmin=0.2, vmax=1)
cbar = ax.figure.colorbar(im, ax=ax, cmap="RdPu")
#plt.clim(0.2, 1)
cbar.ax.set_ylabel('Pearson correlation', rotation=-90, va="bottom")

ax.set_xticks(np.arange(len(timestamps)))
ax.set_yticks(np.arange(len(timestamps)))
ax.set_xticklabels(timestamps)
ax.set_yticklabels(timestamps)
ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

for edge, spine in ax.spines.items():
    spine.set_visible(False)

ax.set_xticks(np.arange(correlation.shape[1]+1)-.5, minor=True)
ax.set_yticks(np.arange(correlation.shape[0]+1)-.5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
ax.tick_params(which="minor", bottom=False, left=False)



if(args.plot):
    plt.savefig(os.path.join(args.plot, "correlation_%s_%s.%s.png" % (args.first[0].split(".")[0],  args.second[0].split(".")[0], args.corrtype)), format = 'png')
else:
    plt.show()

        