#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Checks the recovery and false discovery rate of Peakaboo'''

import argparse
import os
import sys
from itertools import combinations;



import numpy as np;
from scipy.stats import pearsonr, normaltest, variation
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

from afbio.generators import get_only_files



parser = argparse.ArgumentParser(description='Checks the recovery and false discovery rate of Peakaboo');
parser.add_argument('--original', nargs = '?', required=True, type = str, help = "Path to the original peaks, bed file");
parser.add_argument('--detected', nargs = '?', required=True, type = str, help = "Path to the detected annotated peaks, gff format");
parser.add_argument('--maxd', nargs = '?', default=60, type = int, help = "Maximum allowed distance between real and detected peaks");
parser.add_argument('--mode', nargs = '?', required=True, choices = ['length', 'error'], type = str, help = "Type of the variable");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();


def get_name(path, mode):
    if(mode == 'length'):
        return os.path.basename(path).split(".")[0].split("_")[-3]
    if(mode == 'error'):
        return os.path.basename(path).split(".")[0].split("_")[-2]


def find_closest(intervals_1, intervals_2):
    res = []
    start = 0;
    for i1 in intervals_1:
        curint = intervals_2[start]
        curd = abs(i1[0] - curint[0])
        for p, i2 in enumerate(intervals_2[start+1:]):
            d = abs(i1[0] - i2[0])
            if(d<curd):
                curd = d;
                curint = i2
            else:
                res.append((i1, curint, curd))
                start += p
                break;
        else:
            res.append((i1, curint, curd))
            start += p
    return(res)
                
        

def process_detected(path, original, maxd):
    detected = [(int(x.name), float(x.attrs['topcoverage'])) for x in BedTool(path)]
    detected.sort(key = lambda x: x[0])
    recovery = find_closest(detected, original);
    
    true_total = len(original);
    discovered_total = len(detected);
    
    true_positive = len( [ x for x in recovery if x[2]<=maxd] );
    print(true_positive)
    false_positive = discovered_total - true_positive
    
    return true_positive/true_total, 1-false_positive/discovered_total # sensitivity specificity
    
    


original = [(int(x.name), float(x.score)) for x in BedTool(args.original)]
original.sort(key = lambda x: x[0])
#process_detected(args.detected, original, args.maxd)


name2stat = []
for path in [x for x in get_only_files(args.detected) if 'annotated' in x]:
    name = get_name(path, args.mode)
    print(name)
    sens, spec = process_detected(path, original, args.maxd)
    name2stat.append((name, sens*100, spec*100))
    
name2stat.sort(key = lambda x: int(x[0]))




data = [x[1] for x in name2stat],  [x[2] for x in name2stat]
labels = [x[0] for x in name2stat]
fontsize=24 
linewidth = 3 
width = 0.4
adjustmens = [-width/2, width/2]

x = np.arange(len(data[0]))
barcolors = ['lightblue', 'crimson']
barlabels = ['sensitivity', 'specificity']

fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.9])

if(args.mode == 'error'):
    ax.set_xlabel('Error rate [%]', fontsize=fontsize)
if(args.mode == 'length'):
    ax.set_xlabel('Read length [bp]', fontsize=fontsize)
    
ax.set_ylabel('100%', fontsize=fontsize)
ax.set_xticks(x)
ax.set_xticklabels(labels)    
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for axis in ['bottom','left','right']:
    ax.spines[axis].set_linewidth(linewidth)


for adj, yvals, label, color in zip( adjustmens, data, barlabels, barcolors):
    ax.bar(x + adj, yvals, width, label=label, color = color)


fig.legend(loc=(0.25, 0.85), frameon=False, fontsize=fontsize, ncol = 2)
plt.savefig(os.path.join(args.outdir, "sens_spec_%s.%s"  % (args.mode, args.format)) , format = args.format)




#for el in name2stat:
    #print(el)
    
    
    

