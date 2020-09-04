'''Checks co-occurrence of regulators in among multiple bacterial genomes'''

import argparse
import os
import sys
from itertools import combinations
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt









parser = argparse.ArgumentParser(description='Checks co-occurrence of regulators in among multiple bacterial genomes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the occurrence table, csv format");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Path for the output plot");
parser.add_argument('--simnum', nargs = '?', default=100, type = int, help = "Number of random simulations");
#parser.add_argument('--maxlength', nargs = '?', default=10e10, type = int, help = "Maximal length allowed for prophages");
args = parser.parse_args();

def normalized_dot(l1, l2):
    norma = sum(l1)*sum(l2)/len(l1)
    return np.dot(l1, l2)/norma

df = pd.read_csv(args.path, index_col = 0, delimiter = ',')
data = df.to_numpy()
genome_names = list(df.index.values)
regulator_names = list(df.columns.values)

occurrence_list = []
for i in range(data.shape[1]):
    occurrence_list.append([ 1 if x else 0 for x in data[:,i] ])

signal = []
for i, j in combinations(range(len(occurrence_list)), 2):
    signal.append((i, j, normalized_dot(occurrence_list[i], occurrence_list[j])))
                               
noise = []
for _ in range(args.simnum):
    for i, j in combinations(range(data.shape[1]), 2):
        a = occurrence_list[i]
        b = occurrence_list[j]
        noise.append(normalized_dot(random.sample(a, len(a)), random.sample(b, len(b))))

p_levels = [99.99, 99.9, 99, 95]
p_str = [0.0001, 0.001, 0.01, 0.05]
p_vals = [np.percentile(noise, x) for x in p_levels]


cmatrix = np.zeros(( len(occurrence_list), len(occurrence_list) ))
print("\t".join(( 'gene_a', 'gene_b', 'score', 'p_val' )))
for i, j, val in signal:
    cmatrix[i, j] = np.log2(val+1);
    cmatrix[j, i] = np.log2(val+1);
    for p, level in zip(p_vals, p_levels):
        if(val>p):
            print("\t".join((regulator_names[i], regulator_names[j], "%1.1f" % val, "%1.4f" % (1 - level/100) )))
            break;

diagonal = cmatrix.max()
for i in range(cmatrix.shape[0]):
    cmatrix[i,i] = diagonal

        
######################################################################################
### Draw a heatmap

cmap="GnBu"
cmap = 'RdPu'
#cmap="inferno_r"
fontsize = 22
fig, ax = plt.subplots(figsize=(16, 16))
plt.tight_layout(rect=[0.15, 0, 1, 0.8])
im = ax.imshow(cmatrix, cmap=cmap)
cbar = ax.figure.colorbar(im, ax=ax, ticks=[np.log2(x+1) for x in p_vals])
cbar.ax.set_yticklabels(p_str)
cbar.ax.set_ylabel('p-value', rotation=-90, va="bottom", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

ax.set_xticks(np.arange(len(regulator_names)))
ax.set_yticks(np.arange(len(regulator_names)))
ax.set_xticklabels(regulator_names, rotation=90, fontsize=fontsize)
ax.set_yticklabels(regulator_names, fontsize=fontsize)
ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

for edge, spine in ax.spines.items():
    spine.set_visible(False)

ax.set_xticks(np.arange(cmatrix.shape[1]+1)-.5, minor=True)
ax.set_yticks(np.arange(cmatrix.shape[0]+1)-.5, minor=True)
ax.tick_params(which="minor", bottom=False, left=False)

if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format)
else:
    plt.show()




