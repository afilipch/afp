#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Outputs basic statistics on prophages'''

import argparse
import os
import sys
import math
from collections import defaultdict, Counter, namedtuple
from itertools import product, combinations
from glob import glob
import copy


from Bio import Entrez
import numpy as np;

from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;
from matplotlib.patches import Circle
import matplotlib.patches as mpatches



parser = argparse.ArgumentParser(description='Outputs basic statistics on prophages');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the prophage table, tsv format");
parser.add_argument('--maxlength', nargs = '?', default=10e10, type = int, help = "Maximal length allowed for prophages");
parser.add_argument('--btype', nargs = '?', default='Actinobacteria', type = str, help = "Type of bacteria to be analyzed");
parser.add_argument('--maxfraction', nargs = '?', default=0.05, type = float, help = "Maximal allowed prophage genome length (as a fraction of host genome length) to be reported as circle");
args = parser.parse_args();



### Reading the input
prophages = BedTool(args.path)

prophages_filtered = [x for x in prophages if len(x) <= args.maxlength]
sys.stderr.write("%d prophages passed the length filtering (length <= %d) out of %d\n\n" % (len(prophages_filtered), args.maxlength, len(prophages)))
prophages = prophages_filtered;
prophages_filtered = [x for x in prophages if args.btype in x.attrs['host_family']]
sys.stderr.write("%d %s prophages are selected from total %d prophages\n\n" % (len(prophages_filtered), args.btype, len(prophages)))
prophages = prophages_filtered;
prophages_filtered = [x for x in prophages if float(x.attrs['rlength']) < args.maxfraction]
sys.stderr.write("%d prophages with maximum relative length equal %1.3f are selected from total %d prophages\n\n" % (len(prophages_filtered), args.maxfraction, len(prophages)))
prophages = prophages_filtered;


 
scale = 1000;
def prophages2coverage(prophages, scale):
    coverage = np.zeros(scale);
    for prophage in prophages:
        s, e = int(prophage.attrs['host_start']), int(prophage.attrs['host_end'])
        coverage[s:e] += 1;
    return coverage;
        
prophages_intact = [x for x in prophages if x.attrs['integrase'] == 'True']
coverage_intact = prophages2coverage(prophages_intact, scale)


cov_tr = np.mean(coverage_intact)*4
#print(np.mean(coverage_intact))
peak = [];
peaks = [];
for c, x in enumerate(coverage_intact):
    if(x>cov_tr):
        peak.append(x);
    elif(peak):
        #print(c-len(peak), c, max(peak));
        peaks.append((c-len(peak), c));
        peak = [];
else:
    if(peak):
        print(c-len(peak), c, max(peak));
        
        
total = [];
overlapped = [];
for prophage in prophages_intact:
    ps, pe = int(prophage.attrs['host_start']), int(prophage.attrs['host_end']);
    total.append(tuple(prophage.attrs['host_family'].split(",")[2:]));
    for rs, re in peaks:
        if(ps<re and pe>rs):# and tuple(prophage.attrs['host_family'].split(",")[2:]) == ('Corynebacteriales', 'Mycobacteriaceae', 'Mycobacterium', 'Mycobacteriumtuberculosiscomplex')):
            sys.stdout.write(str(prophage))
            overlapped.append(tuple(prophage.attrs['host_family'].split(",")[2:]));
            
            
print()
print(Counter(total)[('Corynebacteriales', 'Mycobacteriaceae', 'Mycobacterium', 'Mycobacteriumtuberculosiscomplex')], len(total));
print(Counter(overlapped)[('Corynebacteriales', 'Mycobacteriaceae', 'Mycobacterium', 'Mycobacteriumtuberculosiscomplex')], len(overlapped));
        
        
    
        
        
        
        
        
        
    





