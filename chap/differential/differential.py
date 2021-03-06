#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Finds consensus regions for the peaks found in different experiments'''

import argparse
import os
import sys
from collections import defaultdict
import copy


import numpy as np;
import pandas as pd;
from pybedtools import BedTool, Interval
from itertools import combinations, permutations
import matplotlib.pyplot as plt;

from afbio.pybedtools_af import construct_gff_interval, intersection2gff
from afbio.peaks import find_shared_peaks, shared_peaks_stat_to_string


parser = argparse.ArgumentParser(description='Finds consensus regions for the peaks found in different experiments');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the merged replicates peaks, gff files");
parser.add_argument('--maxd', nargs = '?', default=60, type = int, help = "Maximal distance allowed between top positions of the peaks");
parser.add_argument('--fraction', nargs = '?', default=1.0, type = float, help = "Minimum fraction of replicates to harbor a valid peak");
parser.add_argument('--mincov', nargs = '?', default=3.0, type = float, help = "Minimum mean (among replicates) coverage for a peak to be considered as expressed greatly than the other");
parser.add_argument('--minfold', nargs = '?', default=2.0, type = float, help = "Minimum fold difference between two peaks to be considered as differential");
parser.add_argument('--area', nargs = '?', default=False, const = True, type = bool, help = "If set, fold change is calculated based on area coverage, ontherwise based on top coverage");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args()


def score_top_coverage(peak):
    topcoverage = [float(x) for x in peak.attrs['topcoverage'].split(",")]
    return np.mean([float(x) for x in topcoverage if x])


def score_area_coverage(peak):
    area_coverage = [float(x) for x in peak.attrs['area_coverage'].split(",")]
    return np.mean([float(x) for x in area_coverage if x])


if(args.area):
    peak_score = score_area_coverage;
else:
    peak_score = score_top_coverage;


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
    
    
    
    


            
######################################################################################################## 

def pair2interval(peak1, peak2):
    diff_peak = copy.copy(peak1);
    if(peak2):
        diff_peak.attrs['other_topcoverage'] = peak2.attrs['topcoverage']
        diff_peak.attrs['other_area_coverage'] = peak2.attrs['area_coverage']
        diff_peak.attrs['other_position'] = peak2.name
        diff_peak.attrs['fold'] = "%1.3f" % (peak_score(peak1)/peak_score(peak2))
    else:
        diff_peak.attrs['other_topcoverage'] = 'None'
        diff_peak.attrs['other_area_coverage'] = 'None'
        diff_peak.attrs['other_position'] = 'None'
        diff_peak.attrs['fold'] = 'None'
    return diff_peak;
        
        
        
    

def check_peak(peak, mincov):
    #print(peak)
    topcoverage = [float(x) for x in peak.attrs['topcoverage'].split(",")]
    return all(topcoverage) and np.mean(topcoverage) >= mincov

def peak_score(peak):
    topcoverage = [float(x) for x in peak.attrs['topcoverage'].split(",")]
    return np.mean([float(x) for x in topcoverage if x])


def compare_two_peaks(peak1, peak2, minfold, mincov):
    if(peak1 and check_peak(peak1, mincov)):
        if(peak2):
            return peak_score(peak1) > minfold*peak_score(peak2);
        else:
            return True 
        
        
def compare_multiple_peaks(compiled, pairs, minfold, mincov):
    temp_d = dict(compiled)
    res = []
    for p in pairs:
        peak1, peak2 = [temp_d.get(x) for x in p]
        if_diff = compare_two_peaks(peak1, peak2, minfold, mincov)
        if(if_diff):
            res.append(( p, peak1, peak2 ))
    return res;
        
        

#def check_replicates(bedtool, fraction_threshold):
    #res = [];
    #for interval in bedtool:
        #topcoverage = [float(x) for x in interval.attrs['topcoverage'].split(",")]
        #fraction = len([x for x in topcoverage if x])/len(topcoverage)
        #if(fraction >= fraction_threshold):
            #res.append(interval);
    #return res;

            
#def find_different(peaks1, peaks2, maxd):
    #diff_peaks = [];
    #for peak in peaks1:
        #close_peaks = [x for x in peaks2 if x.chrom == peak.chrom and abs(int(x.name)-int(peak.name)) <= maxd ]
        #if(not close_peaks):
            #print(peak)
        ##print(len(close_peaks));
    
            

    
        
size = len(args.path)
pairs = list(permutations(range(size), 2))

condition_names = [os.path.basename(x).split(".")[0] for x in args.path]
blist = [BedTool(x) for x in args.path]
res_total, stat_total_counts = find_shared_peaks(blist, args.maxd)
sys.stderr.write(shared_peaks_stat_to_string(stat_total_counts, size))

pos2pairs = defaultdict(list);
for compiled in res_total:
    for p, peak1, peak2 in compare_multiple_peaks(compiled, pairs, args.minfold, args.mincov):
        pos2pairs[p].append((peak1, peak2))
    
for p in pairs:
    label = "%s_GREATER_%s" % (condition_names[p[0]], condition_names[p[1]])
    name = "%s.gff" % label
    sys.stderr.write("%s:\t%d\n" % (label, len(pos2pairs[p])))
    with open(os.path.join(args.outdir, name), 'w') as f:
        for peak1, peak2 in pos2pairs[p]:
            f.write(str(pair2interval(peak1, peak2)))


#condition_bedtools = [check_replicates(BedTool(x), args.fraction) for x in args.path]
    
    
#for (name1, peaks1), (name2, peaks2) in combinations(zip(condition_names, condition_bedtools), 2):
    #find_different(peaks1, peaks2, args.maxd)
    

    
    












