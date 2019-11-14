#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores how motif distance to TSS affects gene expression change'''

import argparse
import os
import sys
#import copy
#from collections import defaultdict, Counter


import numpy as np;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Explores how motif distance to TSS affects gene expression change');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the fimo output file, gff file");
parser.add_argument('--diff', nargs = '?', required=True, type = str, help = "Path to the expression file (0.5h time-point), tsv fomat")
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts regions, gff format");
parser.add_argument('--inside', nargs = '?', default=200, type = int, help = "Maximum allowed distance to TSS while inside a gene");
parser.add_argument('--maxd', nargs = '?', default=800, type = int, help = "Maximum allowed distance to TSS");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Output destination for the plot");
args = parser.parse_args()

def get_center(interval):
    istart, istop = [int(x) for x in interval.name.split(':')[1].split('(')[0].split("-")]
    center = (istart + interval.start) + len(interval)//2
    return center
    #print(interval, center);
    
    


def find_closest(center, starts_plus, stops_minus, gene2diff, d_threshold, inside):
    distances = [('+', x[0], x[1]-center) for x in starts_plus]
    distances.extend([('-', x[0], center-x[1]+1) for x in stops_minus]);
    distances = [x for x in distances if x[2]>-1*inside]
    strand, gene_name, mindistance = min(distances, key = lambda x: abs(x[2]))
    #print(gene_name)
    if(abs(mindistance) <= d_threshold):
        wt, ko, log_ko_wt = gene2diff.get(gene_name, (None, None, 9999999))
        if(log_ko_wt and abs(log_ko_wt)<10 and (ko+wt)>10):
            return mindistance, log_ko_wt, gene_name.split("-")[1]
        else:
            return None
     
        #print(center, gene_name, mindistance);
    
    

gene2diff = {}
with open(args.diff) as f:
    next(f)
    for l in f:
        a = l.strip().split("\t");
        name = a[0]
        wt = np.mean([float(x) for x in a[1].split(";")])
        ko = np.mean([float(x) for x in a[2].split(";")])
        log_ko_wt = float(a[3]);
        gene2diff[name] = (wt, ko, log_ko_wt)
        
        

tr_list = BedTool(args.transcripts)
starts_plus = [ (x.name, x.start) for x in tr_list if x.strand == '+']
stops_minus = [ (x.name, x.stop) for x in tr_list if x.strand == '-']

centers = [get_center(x) for x in BedTool(args.path)]

data = [];
for center in centers:
    data.append(find_closest(center, starts_plus, stops_minus, gene2diff, args.maxd, args.inside))
data = [x for x in data if x]
    
###################### DRAWING SECTION ######################

fontsize=24
xvals = np.array([x[0] for x in data])
yvals = np.array([x[1] for x in data])

fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.95])

ax.set_xlabel('TSS distance', fontsize=fontsize)
ax.set_ylabel('Log2(gene expression fold change 0.5h KO/WT)', fontsize=fontsize)    
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(bottom=-4.5, top=4.5);
ax.scatter(xvals, yvals)
ax.axhline(y=-2, color='red', linewidth=1.5, linestyle='--')
ax.axhline(y=2, color='red', linewidth=1.5, linestyle='--')

prev_y = 0;
ystep = 0.18
for x, y, name in sorted(filter(lambda m: m[1]>2, data), key=lambda m: m[1]):
    pos_x = x+5
    pos_y = y - 0.05
    if(pos_y==prev_y):
        pos_x += 25
    elif(pos_y-prev_y<ystep):
        pos_y =prev_y+ystep
        
    prev_y = pos_y;
    ax.text(pos_x, pos_y, name)


if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format=_format)
else:
    plt.show();


    
