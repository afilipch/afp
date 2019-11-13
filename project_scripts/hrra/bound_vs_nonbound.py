#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Compares genes expression changes for the bound VS unbound targets'''

import argparse
import os
import sys
#import copy
from collections import defaultdict, Counter
from itertools import product;

import numpy as np;
from scipy.stats import sem
#import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Compares genes expression changes for the bound VS unbound targets');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the all-in table");
parser.add_argument('--diff', nargs = 3, required=True, type = str, help = "Path to the file with genes\' differential expression");
parser.add_argument('--distance', nargs = '?', default=400, type = int, help = "Maximum allowed distance for the peak to the closest peak");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path to the plot");
args = parser.parse_args();

mrna_labels = ['0h', '30m', '4h']
#start_dict = {0:0, 1:1, 2:2, 3:2}

def discriminate(gene2diff, bound_set, labels, start):
    bound, unbound = defaultdict(list), defaultdict(list);
    for geneid, currdiff in gene2diff.items():
        for ldiff, label in zip(currdiff[start:], labels[start:]):
            if(sum(ldiff[:2]) > 10):
                change = abs(ldiff[2]);
                if(change < 100):
                    if(geneid in bound_set):
                        bound[label].append(abs(ldiff[2]))
                    else:
                        unbound[label].append(abs(ldiff[2]))
    return bound, unbound;
                    

bound_dict = defaultdict(set);
with open(args.path) as f:
    a = next(f).strip().split("\t")
    chap_labels = ['0h', '30m', '4h']
    #for ck in enumerate(a):
        #print(ck)
    for l in f:
        a = l.strip().split("\t")
        geneid = a[3]
        tss = float(a[6]);
        chap_list = [a[7]] + a[10:12]
        chap_list = [float(x) if x != 'None' else 0 for x in chap_list]
        #FILTERING
        if(tss<=args.distance):
            for chap, label in zip(chap_list, chap_labels):
                if(chap>1):
                    bound_dict[label].add(geneid);
            

gene2diff = defaultdict(list)
for path in args.diff:
    with open(path) as f:
        next(f)
        for l in f:
            a = l.strip().split("\t");
            name = a[0].split("-")[1]
            wt = np.mean([float(x) for x in a[1].split(";")])
            ko = np.mean([float(x) for x in a[2].split(";")])
            log_ko_wt = float(a[3]);
            gene2diff[name].append((wt, ko, log_ko_wt))
       
       
barlabels = []
bound_bars = [];
unbound_bars = [];
bound_errors = [];
unbound_errors = [];

for start, (chap_label, bound_set) in enumerate(bound_dict.items()):
    bound, unbound = discriminate(gene2diff, bound_set, mrna_labels, start)
    for mrna_label in mrna_labels[start:]:
        lb = bound[mrna_label]
        lu = unbound[mrna_label]
        print("%s\t%s\t%1.3f\t%1.3f\t%1.3f\t%1.3f" % (chap_label, mrna_label, np.mean(lb), np.mean(lu), sem(lb), sem(lu)))
        barlabels.append("%s~%s" % (chap_label, mrna_label));
        bound_bars.append(np.mean(lb))
        unbound_bars.append(np.mean(lu))
        bound_errors.append(sem(lb))
        unbound_errors.append(sem(lu))
        
        
        
        
        
############################# DRAWING SECTION #############################
width = 0.4       
fontsize = 24;
linewidth = 4;

fig, ax = plt.subplots(figsize = (16, 9))
x = np.arange(len(bound_bars))
plt.tight_layout(rect=[0.08, 0.08, 1, 1])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

adj = width/2
ax.bar(x, bound_bars, width, label='bound genes', color = 'lightblue')
plt.errorbar(x, bound_bars, yerr=bound_errors, linestyle='', capsize=12, linewidth=linewidth/2, capthick=linewidth/2, color='black')
ax.bar(x+width, unbound_bars, width, label='unbound genes', color = 'coral')
plt.errorbar(x+width, unbound_bars, yerr=unbound_errors, linestyle='', capsize=12, linewidth=linewidth/2, capthick=linewidth/2, color='black')

# add some text for labels, title and axes ticks
ax.set_ylabel('Average absolute Log2 fold change (KO/WT)')
ax.set_xlabel('time of peaks identification ~ time of mRNA expression change')
ax.set_xticks(x+width/2)
ax.legend(frameon=False, fontsize=fontsize, loc ='upper right')
ax.set_xticklabels(barlabels, rotation = 0)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(fontsize)

if(args.plot):
    plt.savefig(args.plot, format = os.path.basename(args.plot).split(".")[-1])
else:
    plt.show()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
