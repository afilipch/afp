#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Draws the chart of expression dynamics over time for phage VS non-phage genes'''

import argparse
import os
import sys
from collections import defaultdict
import copy
import glob


import numpy as np;
from sklearn.cluster import KMeans
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;




parser = argparse.ArgumentParser(description='Draws the chart of expression dynamics over time for phage VS non-phage genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the expression file, tsv format");
parser.add_argument('--minexpr', nargs = '?', default = 50, type = float,  help = "Minimum allowed expression in TPM");
#parser.add_argument('--knum', nargs = '?', default = 6, type = int,  help = "Number of k-means clusters");
parser.add_argument('--size', nargs = '?', default = 50, type = int,  help = "Number of elements per cluster");
parser.add_argument('--outdir', nargs = '?', type = str, help = "Path to the output plot folder");
parser.add_argument('--format', default = 'png', nargs = '?', type = str, help = "Plot format");
args = parser.parse_args()

TIMEPOINTS = 7;
CLUSTER_TYPES = ['before', 'after', 'all']


def transform(datum):
    norma = sum(datum)/100;
    return [x/norma for x in datum]

def split_by_kmeans(data, kmeans):
    res = [[] for _ in range(max(kmeans)+1)]
    for label, datum in zip(kmeans, data):
        res[label].append(datum)
    return res

def get_split_mean(data_split):
    return [np.mean(np.array(x), axis=0) for x in data_split]




transcripts = [];
data = [];
total = 0;
with open(args.path) as f:
    header = next(f).strip().split("\t")
    start = header.index('change') + 1
    xlabels = header[start:start+TIMEPOINTS]
    header.insert(start, 'cluster')

    for l in f:
        total += 1;
        a = l.strip().split("\t")
        #if(a[start-1] == '1'):
        expr = [[float(y) for y in x.split(",")] for x in a[start:start+TIMEPOINTS]]
        #print(expr)
        if( (max([min(x) for x in expr]) >= args.minexpr and all([ max(x) < 2*min(x) + args.minexpr  for x in expr])) or a[start-1] == '1'):
            mean_expr = [np.mean(x) for x in expr]
            if(max(mean_expr[:start]) and max(mean_expr[start:])):
                data.append(mean_expr)
                transcripts.append(a)


#data = np.array(data);
norma = np.sum(np.array(data), axis=0)/1000000;
converted_data = [];
for datum in data:
    converted_data.append([ x[0]/x[1] for x in zip(datum, norma)])
data = converted_data  


data_before = [transform(x[:start]) for x in data]
data_after = [transform(x[start:]) for x in data]
data_all = [transform(x) for x in data]

knum = int(len(data)/args.size) + 1
data_list = [data_before, data_after, data_all]
split_list = []
for data in data_list:
    kmeans = KMeans(n_clusters=knum, random_state=0).fit(data)
    split_list.append(split_by_kmeans(data, kmeans.labels_));

print("\t".join(header));
for tr, c in zip(transcripts, kmeans.labels_):
    tr.insert(start, str(c+1))
    print("\t".join(tr));
        
sys.stderr.write("\ntotal transcripts: %d\ndifferential transcripts: %d\n\n" % (total, len(data)))


if(not args.outdir):
    sys.exit()


for c_type in CLUSTER_TYPES:
    path = os.path.join(args.outdir, c_type)
    if(not os.path.exists(path)):
        os.mkdir(path)
    else:
        for f in glob.glob('%s/*' % path):
            os.remove(f)
        
        
        
#############################################################################################################################
### DRAWING SECTION ###
fontsize = 24
linewidth = 3;
#colors = ('lightblue', 'coral')
#labels = ['non-phage', 'phage']

def draw_trend(data, name, cl_type, xlabels):
    if(cl_type == 'all'):
        xlabs = xlabels
    if(cl_type == 'before'):
        xlabs = xlabels[:start]
    if(cl_type == 'after'):
        xlabs = xlabels[start:]
        
    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])
    for datum in data:
        ax.plot(datum, linewidth=linewidth)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Time", fontsize=fontsize)
    ax.set_ylabel('Expression share [%]', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.set_xticks(range(len(xlabs)))
    ax.set_xticklabels(xlabs)
    ax.text(0.8, 0.8, "size=%d" % len(data), fontsize=fontsize, transform=ax.transAxes)
    
    plt.savefig(os.path.join(args.outdir, cl_type, "%s.%s" % (name, args.format)), format = args.format);
    plt.clf()
    plt.close()
    
 
for data, data_split, cl_type in zip(data_list, split_list, CLUSTER_TYPES):
    draw_trend(data, "all_trends", cl_type, xlabels)
    draw_trend(get_split_mean(data_split), "cluster_trends", cl_type, xlabels)
    for c, cluster_data in enumerate(data_split, start=1):
        draw_trend(cluster_data, "trend_c%d" % c, cl_type, xlabels)

    
