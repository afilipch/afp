#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Selects differentially expressed transcripts for the further experimental investigation'''

import argparse
import os
import sys
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;
import numpy as np
from collections import defaultdict
from scipy.stats import variation



parser = argparse.ArgumentParser(description='Shows an evolution of the peaks for the selected genes over time');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the differential gene expression assigned file, tsv format (differential/sample.assigned.tsv)");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the genomic annotation, plots will be annotated with the provided genomic features")
parser.add_argument('--top', nargs = '?', default=None, type = int, help = "if set, only [--top] candidates are reported");
parser.add_argument('--maxvar', nargs = '?', default=1.01, type = float, help = "Threshold for maximal allowed variation of transcript expression between replicates");
parser.add_argument('--mode', nargs = '?', default='all', choices=['all', 'down', 'up'], type = str, help = "Type of selection, whether all, or only down/up regulated transcripts will be selected");
args = parser.parse_args();


gene2annotation = dict([ (x.name, x.attrs['annotation']) for x in BedTool(args.annotation) ])


def get_varscore(var1, var2, var3, maxvar):
    if(var1 > maxvar or var2 > maxvar):
        return 0;
    else:
        return (var3 - (var1+var2));


def get_score(a, maxvar):
    a1 = [float(x) for x in a[1].split(";")]
    a2 = [float(x) for x in a[2].split(";")]
    var1, var2, var3 = variation(a1), variation(a2), variation([sum(a1), sum(a2)])
    varscore = get_varscore(var1, var2, var3, maxvar)
    change = abs(float(a[3]))
    expression = (sum(a1) + sum(a2))/(len(a1) + len(a2))
    
    return (change**2)*expression*varscore
    
    
    


candidates = []
with open(args.path) as f:
    header = next(f);
    for l in f:
        a = l.strip().split("\t")
        if(args.mode == 'all'):
            candidates.append((a, get_score(a, args.maxvar)))
        elif(args.mode == 'up' and float(a[3])>0):
            candidates.append((a, get_score(a, args.maxvar)))
        elif(args.mode == 'down' and float(a[3])<0):
            candidates.append((a, get_score(a, args.maxvar)))
        else:
            pass;
        
candidates.sort(key=lambda x: x[1], reverse=True)
if(args.top):
    candidates = candidates[:args.top]

header = header.strip() + '\tscore\tannotation' 
print(header)
for cand, score in candidates:
    outline = cand + ["%d" % score, gene2annotation[cand[0]]]
    print("\t".join(outline))
    
    

