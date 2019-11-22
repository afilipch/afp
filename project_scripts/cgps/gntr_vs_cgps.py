#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores how the region of interest changes the chap coverage in varieing conditions'''

import argparse
import os
import sys
import numpy as np;
import pandas as pd;
from collections import defaultdict, Counter
from pybedtools import BedTool

import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Explores how the region of interest changes the chap coverage in varieing conditions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the directory with annotated binding peaks");
parser.add_argument('--interest', nargs = '?', required=True, type = str, help = "Path to the bed file with regions of interest, bed file");
#parser.add_argument('--smooth', nargs = '?', default=0, type = int, help = "Sliding window half-length used to smooth the control genomic coverage, default: 0");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Path to the plot");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output plot directory")
args = parser.parse_args();


peaks_list = [(x.split(".")[0],  BedTool(os.path.join(args.path, x)) ) for x in os.listdir(args.path) if 'annotated' in x];
interest = BedTool(args.interest);

multidata = [];
for pos in range(4):
    data = []
    for name, peaks in peaks_list:
        region = max(peaks.intersect(b=interest, u=True, f=0.35, F=0.35), key = lambda x: float(x.score))
        data.append(( name, float(region.attrs['other_coverage'].split(",")[pos]) ));
    data.sort(key= lambda x: x[0])
    data = [  ( "_".join(x[0][0].split("_")[:2]), x[0][1], x[1][1]) for x in zip(data[::2], data[1::2])]
    data.sort(key= lambda x: x[0][::-1])

    multidata.append(data)

    
    
############################# DRAWING SECTION #############################
def draw(data, ctype, width = 0.6, fontsize = 24, linewidth = 4):
    bars = [(x[1]+x[2])/2 for x in data]
    errors = np.array([(x[1] - min(x[0][1:]), max(x[0][1:]) - x[1]) for x in zip(data, bars)]).transpose()


    fig, ax = plt.subplots(figsize = (16, 9))
    plt.tight_layout(rect=[0.08, 0.08, 1, 1])
    x = np.arange(len(data))
    ax.bar(x, bars, width, label='bound genes', color = 'lightblue')
    plt.errorbar(x, bars, yerr=errors, linestyle='', capsize=12, linewidth=linewidth/2, capthick=linewidth/2, color='black')
    #ax.bar(x+width, [x[2] for x in data], width, label='bound genes', color = 'lightblue')

    # add some text for labels, title and axes ticks
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Name of experiment', fontsize=fontsize)
    ax.set_ylabel('Peak Intensity [averaged coverage]', fontsize=fontsize)
    ax.set_xticks(x)
    ax.set_xticklabels([x[0] for x in data], rotation = 0, fontsize=fontsize-4)
    for item in ax.get_yticklabels():
        item.set_fontsize(fontsize)


    plt.savefig(os.path.join(args.outdir, "interesting_%s.%s") % (ctype, args.format) , format = args.format)

        
        
############################# EXECUTING SECTION #############################
for ctype, data in zip(['all', 'signal', 'noise', 'raw'], multidata):
    draw(data, ctype)
    


