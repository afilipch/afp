#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Draws the chart of expression dynamics over time for phage VS non-phage genes'''

import argparse
import os
import sys
from collections import defaultdict
import copy


import numpy as np;
from sklearn.cluster import KMeans
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;




parser = argparse.ArgumentParser(description='Draws the chart of expression dynamics over time for phage VS non-phage genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the expression file, tsv format");
#parser.add_argument('--minscore', nargs = '?', default = 50000, type = float,  help = "Minimum allowed differential score");
parser.add_argument('--knum', nargs = '?', default = 6, type = int,  help = "Number of k-means clusters");
parser.add_argument('--outdir', required = True, nargs = '?', type = str, help = "Path to the output plot folder");
parser.add_argument('--format', default = 'png', nargs = '?', type = str, help = "Plot format");
args = parser.parse_args()

TIMEPOINTS = 7;
SPLIT_POS = 3;
CLUSTER_TYPES = ['before', 'after', 'all']

for c_type in CLUSTER_TYPES:
    path = os.path.join(args.outdir, c_type)
    if(not os.path.exists(path)):
        os.mkdir(path)

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
    temp = next(f).strip().split("\t")
    start = temp.index('change') + 1
    xlabels = temp[start:start+TIMEPOINTS]
    for l in f:
        total += 1;
        a = l.strip().split("\t")
        if(a[start-1] == '1'):
            data.append([np.mean([float(y) for y in x.split(",")]) for x in a[start:start+TIMEPOINTS]])
            transcripts.append(a[:start])
            #print(a[1])

data_before = [transform(x[:SPLIT_POS]) for x in data]
data_after = [transform(x[SPLIT_POS:]) for x in data]
data_all = [transform(x) for x in data]

#print(data_before)

data_list = [data_before, data_after, data_all]
split_list = []
for data in data_list:
    kmeans = KMeans(n_clusters=args.knum, random_state=0).fit(data)
    split_list.append(split_by_kmeans(data, kmeans.labels_));


        
sys.stderr.write("\ntotal transcripts: %d\ndifferential transcripts: %d\n\n" % (total, len(data)))




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
        xlabs = xlabels[:SPLIT_POS]
    if(cl_type == 'after'):
        xlabs = xlabels[SPLIT_POS:]
        
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

    
