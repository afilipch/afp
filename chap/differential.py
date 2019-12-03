#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Finds consensus regions for the peaks found in different experiments'''

import argparse
import os
import sys
from collections import defaultdict


import numpy as np;
import pandas as pd;
from pybedtools import BedTool, Interval
from itertools import combinations
import matplotlib.pyplot as plt;

from afbio.pybedtools_af import construct_gff_interval, intersection2gff


parser = argparse.ArgumentParser(description='Finds consensus regions for the peaks found in different experiments');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the annotated peaks of the replicates for two different conditions. Conditions must be split by comma, that is 'a_r1 a_r2 a_r3, b_r1 b_r2'");
#parser.add_argument('--rep_numbers', nargs = '+', required=True, type = str, help = "Numbers of replicates per condition. If we provide 3 replicates for condition '1' and the 4 replicates for condition '2' then we have to set '--rep_numbers 3 4'")
parser.add_argument('--zscore', nargs = '?', default=5.0, type = float, help = "Mininimum z-score required for a peak to be a seed for the consensus region");
parser.add_argument('--maxd', nargs = '?', default=60, type = int, help = "Maximal distance allowed between top positions of the peaks");
parser.add_argument('--fraction', nargs = '?', default=1.0, type = float, help = "Minimum fraction of replicates to harbor a valid peak");
#parser.add_argument('--peakflank', nargs = '?', default=5, type = int, help = "Flank's length around a peak to detect maximum coverage");
args = parser.parse_args()


def print_compiled(compiled):
    

def find_shared_double(compiled_dict, current_replicates, maxd):
    next_replicates = [];
    
    r1 = current_replicates[0]
    for r2 in current_replicates[1:]:
        lnext = [];
        for el in r2.intersect(b=r1, wao=True):
            i1, i2 = intersection2gff(el)
            if(i2 and abs(int(i1.name) - int(i2.name)) <= maxd):
                key = len(current_replicates), i2.chrom, int(i2.name)
                if(key not in compiled_dict):
                    compiled_dict[key].append(i2)
                compiled_dict[key].append(i1)
            else:
                lnext.append(i1)
        next_replicates.append(BedTool(lnext));
    #print([len(x) for x in next_replicates])
    return next_replicates
            

    


def find_shared_replicates(replicates, fraction, maxd):
    compiled_dict = defaultdict(list)
    current_replicates = replicates;
    #print()
    
    while(len(current_replicates) > 1):
        #print([len(x) for x in current_replicates])
        current_replicates = find_shared_double(compiled_dict, current_replicates, maxd)
        
    return compiled_dict
        
    #for k, v in compiled_dict.items():
        #if(len(v)==3):
            #print('bu')
            #for el in v:
                #sys.stdout.write(str(el))
            #print()
            #print()
    
        
    


temp = "#".join(args.path).split(",")
temp = [x.strip("#").split("#")  for x in temp]
conditions_list = []
for l in temp:
    conditions_list.append([BedTool(x) for x in l])
    
    
for replicates in conditions_list:
    find_shared_replicates(replicates, args.fraction, args.maxd)
    break
    
    
    












