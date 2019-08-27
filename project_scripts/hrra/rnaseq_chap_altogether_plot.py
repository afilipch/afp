#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Draws a plot for a figure 3'''

import argparse
import os
import sys
from collections import defaultdict, namedtuple
import math

import numpy as np;
import pandas as pd;
from matplotlib import pyplot as plt;
from matplotlib.patches import Rectangle, Arrow
from pybedtools import BedTool, Interval


parser = argparse.ArgumentParser(description='Draws a plot for a figure 3');
parser.add_argument('--regions', nargs = '?', required=True, type = str, help = "Path to the genomic regions of interest, bed format")
parser.add_argument('--chap', nargs = 3, required=True, type = str, help = "Path to the chapseq coverage for time points 0,1,3");
parser.add_argument('--wt', nargs = 3, required=True, type = str, help = "Path to the WT rnaseq coverage for time points 0,1,3");
parser.add_argument('--ko', nargs = 3, required=True, type = str, help = "Path to the KO rnaseq coverage for time points 0,1,3");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the genomic annotation, plots will be annotated with the provided genomic features");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--custom', nargs = '?', default=False, const=True, type = bool, help = "If set the annotation genes are supposed to be already processed, if not they are supposed to be in NCBI gff3 format");
args = parser.parse_args();

def get_cov(path):
    return pd.read_csv(path, sep="\t", names = ["chr", "position", "coverage"]).coverage.values

chap_cov = [get_cov(x) for x in args.chap]
wt_cov = [get_cov(x) for x in args.wt]
ko_cov = [get_cov(x) for x in args.ko]

regions = BedTool(args.regions)
annotation = BedTool(args.annotation)
if(not args.custom):
    annotation = [x for x in annotation if x[2] in ['gene', 'pseudogene']]
#s = 200000
#e = s + 1600
#regions = BedTool([Interval('chr1', s, e, 'test', '0', '+')])

rawannotated = regions.intersect(annotation, wo = True)
region2annotation = defaultdict(list);
for el in rawannotated:
    an = max(el.start, int(el[9])),  min(el.end, int(el[10])),  dict( [x.strip().split('=') for x in el[14].split(";")])['Name'], el[12]
    region2annotation[el.name].append(an)
    
#wt = [x[s:e] for x in wt_cov]
#ko = [x[s:e] for x in ko_cov]
#chap = [x[s:e] for x in ko_cov]
#locan = region2annotation[regions[0].name]    


def draw_annotation(ax, locan, color, start, end):
    for c, an in enumerate(locan):
        rect = Rectangle( (an[0], 0.25+c), an[1] - an[0], 0.5, facecolor = color, edgecolor = color)
        l = min((end-start)/20, an[1]-an[0]);
        if(an[3] == '+'):
            arrow = Arrow(an[1]-l, 0.5+c, l, 0, width=.75, facecolor = 'black', edgecolor = 'black')
        else:
            arrow = Arrow(an[0]+l, 0.5+c, -l, 0, width=.75, facecolor = 'black', edgecolor = 'black')
        ax.add_patch(rect)
        ax.add_patch(arrow)


def setticks(start, end):
    l = float(end - start);
    a = [1, 2, 5,10,20,50,100,200,500,1000, 2000]
    for el in a:
        if( 3 < l/el < 9):
            scale = el;
            break;
    else:
        return [];
    
    
    s = start//scale*scale
    if(s < start):
        s += scale;
        
    e = end//scale*scale
    if(e == end):
        e += scale;
    
        
    locs = [x for x in range(s, e+scale, scale)]
    labels = [str(x) for x in locs]
    return locs, labels


def doublebar(ax, d1, d2, prange, ylim, start, end):
    #print(type(ax))
    ax.bar(prange, d1, 1, color='lightblue');
    ax.bar(prange, -d2, 1, color='darkblue');
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([]);
    ax.set_ylim(*ylim)
    ax.set_xlim(start, end)
    ax.tick_params(axis='both', labelsize='xx-large')
    
def singlebar(ax, d, prange, ylim, start, end):
    ax.bar(prange, d, 1, color=(242/256, 97/256, 68/256));
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([]);
    ax.set_ylim(*ylim)
    ax.set_xlim(start, end)
    ax.tick_params(axis='both', labelsize='xx-large')
        


def wholeplot(chap, rseq1, rseq2, locan, start, end, output):
    prange = range(start, end);
    size_rseq = len(rseq1)
    size_chap = len(chap)
    size = size_chap + size_rseq
    size_an = len(locan)
    fig, axes = plt.subplots(nrows=size, ncols = 1, sharex=False, sharey=False, figsize = (16, 5*(size+1)), frameon=False)
    #plt.suptitle("Top %s peak with z-score %s" % (ordinal(fignum), zscore), fontsize = 'xx-large')
    plt.tight_layout(rect=[0.05, 0.02, 0.95, 0.92 - 0.02*size_an], h_pad = 2)
    
    ###CHAP seq 
    ylim = 0, max([max(x) for x in chap])
    for ax, d in zip(axes, chap): 
        singlebar(ax, d, prange, ylim, start, end);
    
    ###RNA seq
    ylim = max([max(x) for x in rseq2] + [max(x) for x in rseq1])
    ylim = -ylim, ylim
    for ax, d1, d2 in zip(axes[size_chap:], rseq1, rseq2): 
        doublebar(ax, d1, d2, prange, ylim, start, end);
        
    box = axes[0].get_position()
    x0 = box.xmin
    x1 = box.xmax-x0
    ylen = box.ymax - box.ymin
    anax = fig.add_axes([x0, 0.94 - 0.02*size_an, x1, 0.02*size_an])
    anax.set_xlim(start, end)
    
    anax.spines['bottom'].set_visible(False)
    anax.spines['right'].set_visible(False)
    anax.spines['left'].set_visible(False) 
    
    anax.tick_params(axis='both', labelsize='xx-large')
    anax.set_ylim(0, size_an)
    
    xticklocs, xticklabels = setticks(start, end)
    anax.set_xticks(xticklocs)
    anax.set_xticklabels(xticklabels)
    anax.set_yticks([0.5 + x for x in range(size_an)])
    anax.set_yticklabels([x[2] for x in locan], fontsize='xx-large')
    anax.xaxis.tick_top()
    anax.set_xlabel('base position', fontsize='xx-large')    
    anax.xaxis.set_label_position('top') 
    
    draw_annotation(anax, locan, '0.75', start, end);#max([max(x[0]) for x in signal_noise_local]) )
    anax.set_xlim(start, end)
    
    plt.savefig(output, format = args.format)
    
    

for region in regions:
    start, end = region.start, region.end
    print(start, end)
    wt = [x[start:end] for x in wt_cov]
    ko = [x[start:end] for x in ko_cov]
    chap = [x[start:end] for x in chap_cov]
    locan = region2annotation[region.name]   
    wholeplot(chap, wt, ko, locan, start, end, os.path.join(args.outdir, "%s.%s" % (region.name, args.format)));
    
    
