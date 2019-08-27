#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Assigns differential expression when at least one of two experiments is performed  without replicates'''

import argparse
import os
import sys
from collections import defaultdict
import math

import numpy as np;
from scipy.stats import variation 
import matplotlib.pyplot as plt;
from scipy.stats.stats import pearsonr;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Assigns differential expression when at least one of two experiments is performed  without replicates');
parser.add_argument('--first', nargs = '+', required=True, type = str, help = "Path to gene expression data (replicates) for the first sample");
parser.add_argument('--second', nargs = '+', required=True, type = str, help = "Path to gene expression data (replicates) for the second sample");
parser.add_argument('--labels', nargs = 2, required=True, type = str, help = "Names of the samples");
parser.add_argument('--minexpr', nargs = '?', default=0.001, type = float, help = "Minimum required expression for the differential analysis");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();


#def get_sampleid(multipath):
    #return os.path.basename(multipath[0]).split(".")[0]

def get_expr(mp1, mp2):
    genes2expr = defaultdict(lambda: defaultdict(list));
    for path in mp1:
        for interval in BedTool(path):
            genes2expr[interval.name][0].append(float(interval.attrs['tpm']))
    for path in mp2:
        for interval in BedTool(path):
            genes2expr[interval.name][1].append(float(interval.attrs['tpm']))
    return genes2expr;


genes2expr = [ (x[0], x[1][0], x[1][1], sum(x[1][0]), sum(x[1][1]) ) for x in get_expr(args.first, args.second).items()]
genes2expr = [x for x in genes2expr if x[3] > args.minexpr and x[4] > args.minexpr]
genes2expr.sort(key= lambda x: x[3], reverse=True)

#print(genes2expr)
header = "gene_id\t%s\t%s\tlog2(%s/%s)\tnormed_fold" % (args.labels[0], args.labels[1], args.labels[0], args.labels[1]);
print(header)
for gene, tpm1, tpm2, e1, e2 in genes2expr:
    print('%s\t%s\t%s\t%1.3f\t%1.3f' % (gene, ";".join([str(x) for x in tpm1]), ";".join([str(x) for x in tpm2]), math.log2(e1/e2), (e1-e2)/(e1+e2)  ))
    
    
    
###generate loglog plot
scatter = np.array([ (math.log2(x[3]+x[4]), math.log2(x[3]/x[4])) for x in genes2expr]);


fig, ax = plt.subplots(figsize=(16, 9))
plt.plot(scatter[:,0], scatter[:,1], marker = 'o', markerfacecolor='lightblue', markeredgecolor='lightblue', linestyle = "None", markersize = 4);

plt.ylabel("log2(%s/%s)" % tuple(args.labels));
plt.xlabel("log2(TPM)")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#plt.title("Expression vs Fold Change")

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize('xx-large')


plt.savefig(os.path.join(args.outdir, '%s_%s.scatter.%s' % (args.labels[0], args.labels[1], args.format)), format = args.format)
