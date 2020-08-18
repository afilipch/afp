import argparse
import os
import sys
import numpy as np;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

from afbio.pybedtools_af import read_comments




parser = argparse.ArgumentParser(description='Compares expression for the genes coming from different origins dependent on experimental conditions');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the binding peaks, bed/gff format. ORDER: 'genomic', 'cgp3', 'cl31'");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
args = parser.parse_args();

TTYPES = ['genomic', 'cgpN', 'cl31']
barlabels = read_comments(args.path[0])[0].split("=")[1].split(",")
data = []

for path in args.path:
    temp = [];  
    for tr in BedTool(path):
        expr = [np.mean([float(y) for y in x.split(",")]) for x in tr.attrs['expression'].split(":")]
        temp.append(expr)
    l = np.array(temp).sum(axis=0)
    data.append(l)
    
data = np.array(data)
data = data/data.sum(axis=0,keepdims=1)
#data = data.transpose()


###Plotting all data in one plot
    
fontsize=28
linewidth = 5 
width = 0.2
adjustmens = [-width*1.5, -width/2, width/2, width*1.5]

x_range = np.arange(len(data[0]))
barcolors = ['coral', 'lightblue', 'limegreen']

fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.9])

ax.set_xlabel('origin of transcripts', fontsize=fontsize)
ax.set_ylabel('share of expression', fontsize=fontsize)
ax.set_xticks(x_range)
ax.set_xticklabels(barlabels)    
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for axis in ['bottom','left']:
    ax.spines[axis].set_linewidth(linewidth)

curbottom = np.zeros(len(data[0]))
for yvals, label, color in zip(data, TTYPES, barcolors):
    plt.bar(x_range, yvals, width, color=color, label=label, bottom = curbottom)
    curbottom += yvals



fig.legend(loc=(0.25, 0.85), frameon=False, fontsize=fontsize, ncol = 2)
plt.savefig(os.path.join(args.outdir, "condition_by_origin.%s"  % args.format) , format = args.format)

    




 
    






    






























