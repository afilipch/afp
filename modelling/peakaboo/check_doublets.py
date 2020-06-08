#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Checks the recovery and false discovery rate of Peakaboo'''

import argparse
import os
import sys
from itertools import combinations;
from collections import defaultdict
import copy


import numpy as np;
from scipy.stats import pearsonr, normaltest, variation
import matplotlib.pyplot as plt;
from pybedtools import BedTool

from afbio.sequencetools import split2chunks
from afbio.generators import get_only_files


parser = argparse.ArgumentParser(description='Checks the recovery and false discovery rate of Peakaboo');
parser.add_argument('--original', nargs = '?', required=True, type = str, help = "Path to the original peaks, bed file");
parser.add_argument('--detected', nargs = '?', required=True, type = str, help = "Path to the directory with detected annotated peaks, gff format");
parser.add_argument('--maxd', nargs = '?', default=60, type = int, help = "Maximum allowed distance between real and detected peaks");
parser.add_argument('--detection_fraction', nargs = '?', default=75, type = float, help = "Minimum allowed fraction of the detected double peaks");
parser.add_argument('--mode', nargs = '?', required=True, choices = ['length', 'ratio'], type = str, help = "Type of the variable");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();


def get_name(path, mode):
    if(mode == 'length'):
        return os.path.basename(path).split(".")[0].split("_")[-3]
    if(mode == 'ratio'):
        return os.path.basename(path).split(".")[0].split("_")[0][8:]
    

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
                
        

def process_detected(path, original):
    detected = [(int(x.name), float(x.attrs['topcoverage'])) for x in BedTool(path)]
    detected.sort(key = lambda x: x[0])
    recovery = find_closest(original, detected);
    
    length2recovery = defaultdict(list)
    for r1, r2 in split2chunks(recovery, 2):
        if(r1[2]<=args.maxd and r2[2]<=args.maxd):
            length2recovery[r2[0][0] - r1[0][0]].append( int(r1[1] != r2[1]) )
        elif(r1[2]<=args.maxd or r2[2]<=args.maxd):
            length2recovery[r2[0][0] - r1[0][0]].append(0)
            

    xticks, yvals = [], []
    for k, v in sorted(length2recovery.items(), key = lambda x: x[0]):
        xticks.append(k)
        yvals.append( sum(v)/len(v)*100)
        
    return xticks, yvals
   
   
   
def find_step_for_sample(yvals, xticks, name, fraction):
    for y, x in zip(yvals, xticks):
        if(y > fraction):
            return x
    else:
        return int(name)
   
   
### DRAWING SECTION ###

def draw_single(yvals, xticks, name):
    xvals = range(len(yvals))

    fontsize = 20
    fig, ax = plt.subplots(figsize=(16,9))
    ax.set_xlabel('Double peak distance [bp]', fontsize=fontsize)
    ax.set_ylabel("Accuracy [%]", fontsize=fontsize)    
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.bar(xvals, yvals, 0.98, color = 'darkblue')   
    plt.xticks(xvals, xticks, rotation = 45)

    plt.savefig(os.path.join(args.outdir, "%s.%s"  %  (name, args.format) ) , format = args.format)
    plt.clf()
    plt.close()




if(args.mode == 'length'):
    original = [(int(x.name), float(x.score)) for x in BedTool(args.original)]
    original.sort(key = lambda x: x[0])
    #for r1, r2 in split2chunks(original, 2):
        #print(abs(r1[0] - r2[0]))

    length2step = []
    for path in [x for x in get_only_files(args.detected) if 'annotated' in x]:
        name = get_name(path, args.mode)
        xticks, yvals = process_detected(path, original);
        length2step.append(( int(name),  find_step_for_sample(yvals, xticks, name, args.detection_fraction) ))
        draw_single(yvals, xticks, name);
    xlabel = 'Read length'  
   
   
elif(args.mode == 'ratio'):
    original_dict = {}
    for path in get_only_files(args.original):
        name = os.path.basename(path).split(".")[0][1:]
        original = [(int(x.name), float(x.score)) for x in BedTool(path)]
        original.sort(key = lambda x: x[0])
        original_dict[name] = copy.copy(original);
        
    length2step = []
    for path in [x for x in get_only_files(args.detected) if 'annotated' in x]:
        name = get_name(path, args.mode)
        original = original_dict[name]
        xticks, yvals = process_detected(path, original);
        length2step.append(( int(name),  find_step_for_sample(yvals, xticks, name, args.detection_fraction) ))
        draw_single(yvals, xticks, name);
        if(name == '4'):
            detected = [(int(x.name), float(x.attrs['topcoverage'])) for x in BedTool(path)]
            detected.sort(key = lambda x: x[0])
            recovery = find_closest(original, detected);
            #print(original)
            
            length2recovery = defaultdict(list)
            for r1, r2 in split2chunks(recovery, 2):
                #print(r1);
                #print(r2)
                #print()
                if(r1[2]<=args.maxd and r2[2]<=args.maxd):
                    length2recovery[r2[0][0] - r1[0][0]].append( int(r1[1] != r2[1]) )
                elif(r1[2]<=args.maxd or r2[2]<=args.maxd):
                    length2recovery[r2[0][0] - r1[0][0]].append(0)
            print(length2recovery)
        
    xlabel = 'Double peak ratio'     
        
        
        
### Drawing All together    
length2step.sort(key = lambda x: x[0])  
xvals = [x[0] for x in length2step]
yvals = [x[1] for x in length2step]
    
fontsize = 20
fig, ax = plt.subplots(figsize=(16,9))
ax.set_xlabel(xlabel, fontsize=fontsize)
ax.set_ylabel("Double peak distance [bp]", fontsize=fontsize)    
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
plt.xticks(xvals, ['%d' % x for x in xvals])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.plot(xvals, yvals, color = 'darkblue', marker='o')   
plt.savefig(os.path.join(args.outdir, "all.%s"  %   args.format ) , format = args.format)


