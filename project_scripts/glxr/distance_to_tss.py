#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates table with all the information regarding GLXR camp project'''

import argparse
import os
import sys

import numpy as np;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;




parser = argparse.ArgumentParser(description='Draws a TSS distribution plot');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated consensus regions");
parser.add_argument('--mincov', nargs = '?', default=1, type = float, help = "Minimum allowed coverage")
parser.add_argument('--format', nargs = '?', default='svg', type = str, help = "Plot format, svg by default");
parser.add_argument('--outdir', nargs = '?', required=False, type = str, help = "Path to the output plot directory")
args = parser.parse_args();


def fix_density(ax, bins, scale):
    factor = bins[1]-bins[0]
    ypos = [float(x) for x in ax.get_yticks().tolist()]
    yvals = [float(x)*factor for x in ax.get_yticks().tolist()]
    
    step = scale/yvals[1]*ypos[1]
    newvals = [x for x in np.arange(0, max(yvals), scale)]
    newpos = np.arange(0, step*len(newvals), step)

    #print(sum(ylabels))
    ax.set_yticks(newpos)
    ax.set_yticklabels(["%d" % (x*100) for x in newvals])


def check(region, mincov):
    return len([x for x in region.attrs['topcoverage'].split(",")[:3] if float(x)>mincov ]) > 1;

annpeaks = [x for x in BedTool(args.path) if check(x, args.mincov)]


fontsize=24
linewidth = 5
scores = [float(x.attrs['tss']) for x in annpeaks if x.attrs['tss'] != 'nan']
scores.sort();
selected_scores = [x for x in scores if x <= 300 and x >= -100]
#print(min(scores))

fig, axes = plt.subplots(ncols=2, figsize = (22, 7), frameon=False)
fig.tight_layout(rect=[0.05, 0.1, 1, 1])
fig.subplots_adjust(wspace = 0.2)
for data, ax in zip([scores, selected_scores], axes):
    _, bins, _ = ax.hist(data, bins = 20, density = True)
    
    ax.set_xlabel('TSS distance', fontsize=fontsize)
    ax.set_ylabel('Fraction [%]', fontsize=fontsize)    
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fix_density(ax, bins, 0.05)
    

plt.savefig(os.path.join(args.outdir, "tss_hist_selected.%s") % args.format, format = args.format)
plt.close()
