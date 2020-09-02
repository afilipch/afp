#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores various aspects of interplay between binding motifs and binding peaks'''

import argparse
import os
import sys
from collections import defaultdict
import copy
from pathlib import Path


import numpy as np;
from pybedtools import BedTool
from itertools import combinations, permutations
import matplotlib.pyplot as plt;
from afbio.pybedtools_af import read_comments
from afbio.peaks import find_closest_feature_unstranded, find_proximal_feature_unstranded, find_closest_peak_unstranded




parser = argparse.ArgumentParser(description='Explores various aspects of interplay between binding motifs and binding peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated (merged) peaks, gff format");
parser.add_argument('--fimo', nargs = '?', required = True, type = str, help = "Path to the fimo output, gff format");
parser.add_argument('--maxd', nargs = '?', default = 50, type = int, help = "Minimal distance from motif to the closest peaks");
parser.add_argument('--outdir', nargs = '?', required = True, type = str, help = "Path to the output directory");
parser.add_argument('--format', nargs = '?', default = 'png', type = str, help = "Plots format");
parser.add_argument('--color', nargs = '?', default = 'lightblue', type = str, help = "Color of the plots");
args = parser.parse_args()

def zscore_vs_motif(peaks, labels, score_type):
    mylabels = ['mean'] + labels
    result = [[] for _ in range(len(mylabels))]
    for peak in peaks:
        ifo = min(float(peak.attrs['motif']),  1)
        zscores = [float(x) for x in peak.attrs[score_type].split(',')]
        result[0].append((np.mean(zscores), ifo))
        for p, zs in enumerate(zscores, start = 1):
            if(zs):
                result[p].append((zs, ifo))
    for mylist in result:
        mylist.sort(key = lambda x: x[0], reverse = True)
        
    return mylabels, result


def draw_zscore_vs_motif(data, path, format_, color, limits, xlabel, fontsize=24, linewidth=3):
    yvalues = []
    counts = []
    for l1, l2 in zip(limits, limits[1:]):
        lf = [ x[1] for x in data if x[0]>l1 and x[0]<=l2]
        if(lf):
            yvalues.append(np.mean(lf))
        else:
            yvalues.append(0)
        counts.append(len(lf))
    lf = [ x[1] for x in data if x[0]>limits[-1]]
    if(lf):
        yvalues.append(np.mean(lf))
    else:
        yvalues.append(0)
    counts.append(len(lf))
    
    xlabels = ["%d~%d" % x for x in zip(limits, limits[1:])] + [">%d" % limits[-1]]
            
    fig, ax = plt.subplots(figsize=(16,9))
    #plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Peak fraction with motifs', fontsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    
    brange = range(len(yvalues))
    ax.bar(brange, yvalues, 0.75, color=color)
    plt.xticks(brange, xlabels)
    for x, y, c in zip(brange, yvalues, counts):
        ax.text(x, y, str(c), horizontalalignment='center', verticalalignment='bottom', fontsize=fontsize)
    plt.savefig(path, format = format_);
    plt.close()
    
    
def replicates_vs_motif(peaks, labels):
    result = [[] for _ in range(len(labels))]
    for peak in peaks:
        ifo = min(float(peak.attrs['motif']),  1)
        zscores = [float(x) for x in peak.attrs['zscores'].split(',')]
        replicates = len([x for x in zscores if x])
        result[replicates-1].append(ifo)
    return result
        
    
    
def draw_replicates_vs_motif(data, path, format_, color, fontsize=24, linewidth=3):
    yvalues = [np.mean(x) for x in data]
    counts = [len(x) for x in data]
    brange = range(len(yvalues))
    xlabels = [str(x+1) for x in brange]
    
            
    fig, ax = plt.subplots(figsize=(16,9))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Peak fraction with motifs', fontsize=fontsize)
    ax.set_xlabel('Number of replicates', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    
    ax.bar(brange, yvalues, 0.75, color=color)
    plt.xticks(brange, xlabels)
    for x, y, c in zip(brange, yvalues, counts):
        ax.text(x, y, str(c), horizontalalignment='center', verticalalignment='bottom', fontsize=fontsize)    
    
    plt.savefig(path, format = format_);
  
  
  
def num_motifs(peaks):
    result = defaultdict(int);
    for peak in peaks:
        ifo = int(peak.attrs['motif'])
        result[ifo] += 1
    return [x[1] for x in sorted(result.items(), key = lambda x: x[0])]
    
    
def draw_num_motifs(yvalues, path, format_, color, fontsize=24, linewidth=3):
    brange = range(len(yvalues))
    xlabels = [str(x) for x in brange]
            
    fig, ax = plt.subplots(figsize=(16,9))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Number of motifs', fontsize=fontsize)
    ax.set_xlabel('Number of motifs per peak', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    
    ax.bar(brange, yvalues, 0.75, color=color)
    plt.xticks(brange, xlabels)
    
    plt.savefig(path, format = format_);
    
    
def draw_closest_motifs(data, path, format_, color, fontsize=24, linewidth=3):
    fig, ax = plt.subplots(figsize=(16,9))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Number of motifs', fontsize=fontsize)
    ax.set_xlabel('Number of motifs per peak', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    
    ax.hist(data, bins = np.linspace(-60, 60, 13), color=color)
    
    plt.savefig(path, format = format_);
    

zscore_fraction_dir = os.path.join(args.outdir, 'zscore_fraction')
Path(zscore_fraction_dir).mkdir(parents=True, exist_ok=True)
topcoverage_fraction_dir = os.path.join(args.outdir, 'topcoverage_fraction')
Path(topcoverage_fraction_dir).mkdir(parents=True, exist_ok=True)

comments = read_comments(args.path)[0]
labels = comments.split(' ')[1].split(",");
fimo = BedTool(args.fimo)
peaks = BedTool(args.path)
newpeaks = []

with open(os.path.join(args.outdir, 'annotated.gff') , 'w') as f:
    f.write("%s\n" % comments)
    for peak, num in find_proximal_feature_unstranded(peaks, fimo, args.maxd):
        peak.attrs['motif'] = str(num)
        newpeaks.append(peak)
        f.write(str(peak))

            
with open(os.path.join(args.outdir, 'orphans.gff'), 'w') as f:
    for feature, peak, distance in find_closest_peak_unstranded(peaks, fimo):
        if(abs(distance)>args.maxd*2):
            f.write(str(feature));


##################################################################################################################################
### Drawing section ###

score_types = ('zscores', 'topcoverage')
limits_list = ([0, 1, 2, 5, 20, 100], [0, 1, 2, 3, 5, 10, 20])
xlabels = ['Z-score', 'Top coverage']
dirs = (zscore_fraction_dir, topcoverage_fraction_dir)

for score_type, limits, xlabel, curdir in zip(score_types, limits_list, xlabels, dirs):
    f_labels, f_data = zscore_vs_motif(newpeaks, labels, score_type)
    for label, data in zip(f_labels, f_data):
        path = os.path.join(curdir, '%s.%s' % (label, args.format))
        draw_zscore_vs_motif(data, path, args.format, args.color, limits, xlabel, fontsize=24, linewidth=3)
    
r_data = replicates_vs_motif(newpeaks, labels)
draw_replicates_vs_motif(r_data, os.path.join(args.outdir, "replicates_vs_motif.%s" % args.format) , args.format, args.color, fontsize=24, linewidth=3)

n_data = num_motifs(newpeaks)
draw_num_motifs(n_data, os.path.join(args.outdir, "num_motifs.%s" % args.format) , args.format, args.color, fontsize=24, linewidth=3)

distances = [x[2] for x in find_closest_feature_unstranded(peaks, fimo)]
draw_closest_motifs(distances, os.path.join(args.outdir, "distance_motifs.%s" % args.format) , args.format, args.color, fontsize=24, linewidth=3)

    


