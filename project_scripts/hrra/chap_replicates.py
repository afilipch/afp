#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Augments the super table with the information on how frequently a peak appears in the replicates'''

import argparse
import os
import sys

import numpy as np;
from pybedtools import BedTool, Interval
from afbio.generators import get_only_files
from collections import defaultdict


parser = argparse.ArgumentParser(description='Augments the super table with the information on how frequently a peak appears in the replicates');
parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Path to the super table with peaks")
parser.add_argument('--replicates', nargs = '+', required=True, type = str, help = "Path to the folders with binding peaks");
parser.add_argument('--journal', nargs = '?', required=True, type = str, help = "Path to the super table in journal format")
args = parser.parse_args();


name2string = {};
with open(args.journal) as f:
    next(f);
    for l in f:
        a = l.strip().split("\t")
        name =  "_".join([a[4], a[5]] + a[7:12])
        name2string[name] = a;
        #print(name)
        
#print(len(name2string))
        

        



reference = [];
with open(args.table) as f:
    next(f);
    for l in f:
        a = l.strip().split("\t")
        start = int(a[1])
        stop = int(a[2])
        temp = ["%1.1f" % float(x) if x!='None' else '0.0' for x in a[9:14]]
        interval = Interval("NC_003450.3", start, stop, "_".join(a[5:7] + temp), '0', '+' )
        #print(interval.name)
        reference.append(interval)
        
reference = BedTool(reference);

overlap_counter = defaultdict(int)
for folder in args.replicates:
    files = [x for x in get_only_files(folder) if "filtered" in x]
    replicates_list = [BedTool(x) for x in files];
    for r in reference.intersect(b=replicates_list, u = True):
        overlap_counter[(r.name)] += 1;
    #for km in replicates_list:
        #for r in reference.intersect(b=km, u = True):
            #overlap_counter[(r.start, r.stop)] += 1;
 
for name, a in name2string.items():
    #if(overlap_counter[name]==2):
    a.append(str(overlap_counter[name]))
    print("\t".join(a))
    
    
#for k, v in overlap_counter.items():
    #m = name2string[k]
    
        
#temp = list(overlap_counter.items())
#temp.sort(key = lambda x: x[0])
#for el in temp:
    #print(el);
    
    
    
