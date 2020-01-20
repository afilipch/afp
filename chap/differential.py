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
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the merged replicates peaks, gff files");
parser.add_argument('--maxd', nargs = '?', default=60, type = int, help = "Maximal distance allowed between top positions of the peaks");
parser.add_argument('--fraction', nargs = '?', default=1.0, type = float, help = "Minimum fraction of replicates to harbor a valid peak");
args = parser.parse_args()



def print_compiled(compiled, size):
    temp_d = dict(compiled)
    compiled_processed = [temp_d.get(x, None) for x in range(size)]
    
    area_coverage = ",".join([x.attrs['area_coverage'] if x else '0' for x in compiled_processed])
    topcoverage=",".join([x.attrs['topcoverage'] if x else '0' for x in compiled_processed])
    compiled_processed = [x for x in compiled_processed if x]
    
    compiled = [x[1] for x in sorted(compiled, key = lambda x: x[0])]
    pos = int(sum([int(x.name)*float(x.attrs['topcoverage']) for x in compiled_processed])/sum([float(x.attrs['topcoverage']) for x in compiled_processed]))
    start = min([x.start for x in compiled_processed])
    stop = min([x.stop for x in compiled_processed])
    
    consensus = construct_gff_interval(compiled[0].chrom, start, stop, 'consensus', score='0', strand='.', source='.', frame='.', attrs=[('Name', pos), ('topcoverage', topcoverage), ('area_coverage', area_coverage)])
    
    return consensus;
    
    #print(consensus)
    #for cs in compiled_processed:
        #print(cs)
    #print()
    #print("*"*150)
    
    
    


            
########################################################################################################

def check_replicates(bedtool, fraction_threshold):
    res = [];
    for interval in bedtool:
        topcoverage = [float(x) for x in interval.attrs['topcoverage'].split(",")]
        fraction = len([x for x in topcoverage if x])/len(topcoverage)
        if(fraction >= fraction_threshold):
            res.append(interval);
    return res;

            
def find_different(peaks1, peaks2, maxd):
    diff_peaks = [];
    for peak in peaks1:
        close_peaks = [x for x in peaks2 if x.chrom == peak.chrom and abs(int(x.name)-int(peak.name)) <= maxd ]
        if(not close_peaks):
            print(peak)
        #print(len(close_peaks));
    
            

    
        
    
condition_names = [os.path.basename(x).split(".")[0] for x in args.path]
condition_bedtools = [check_replicates(BedTool(x), args.fraction) for x in args.path]
    
    
for (name1, peaks1), (name2, peaks2) in combinations(zip(condition_names, condition_bedtools), 2):
    find_different(peaks1, peaks2, args.maxd)
    

    
    












