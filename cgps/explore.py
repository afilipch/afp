#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores various aspects of cgps binding'''

import argparse
import os
import sys
#import scipy
import numpy as np;
import pandas as pd;
from pybedtools import BedTool, Interval
from Bio import SeqIO
import heapq

from afbio.sequencetools import sliding_window

#from afbio.filters import dsk, usk, trk;
#from afbio.peaks import convolute


parser = argparse.ArgumentParser(description='Explores various aspects of cgps binding');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks");
parser.add_argument('--coverage', nargs = '?', required=True, type = str, help = "Path to the binding coverage track");
parser.add_argument('--gctrack', nargs = '?', required=True, type = str, help = "Path to the GC track");
parser.add_argument('--gcdrops', nargs = '?', required=True, type = str, help = "Path to the GC drops");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format");


parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");

args = parser.parse_args();

FLANK = 50;

regions = BedTool(args.path);
gc_content = pd.read_csv(args.gctrack, sep="\t" , names = ["chr", "position", "gc"]).gc.values
coverage = pd.read_csv(args.coverage, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values

#read gc extrema
minima, maxima = [0]*len(coverage), [0]*len(coverage)
gcdrops = [];
with open(args.gcdrops) as f:
    current_drop = []
    for l in f:
        a = l.strip().split("\t");
        if(not current_drop and a[3] == 'max'):
            current_drop.append(int(a[1]));
        if(current_drop and a[3] == 'min'):
            current_drop.append(int(a[1]));
            gcdrops.append(tuple(current_drop));
            current_drop = [];
        #if(a[3] == 'min'):
            #minima[int(a[1])] = float(a[2]);
        #else:
            #maxima[int(a[1])] = float(a[2]);
            
            
gcdrops_intervals = []             
for c, (start, end) in enumerate(gcdrops, start=1):
    score = min(gc_content[start: end])
    if(score <= 0.35):
        gcdrops_intervals.append(Interval('chr1', start, end, name = "drop%d" % c, score = "%1.5f"  % score, strand = '+'));

gcdrops_intervals = BedTool(gcdrops_intervals);        
print(len(gcdrops_intervals))

gcdrops_regions = gcdrops_intervals.intersect(regions, c = True)
#for interval in gcdrops_regions:
    #print(interval)


        
            
            
#gcmeans = [np.mean(x) for x in sliding_window(gc_content, 20)]
#print(heapq.nsmallest(20, gcmeans))

#print(heapq.nsmallest(20, minima))  

#selected = [x for x in regions if float(x.score) > 450]
#for region in selected:
    #lmin = min(minima[region.start-FLANK:region.end+FLANK])
    #lmin2 = min(minima[region.start:region.end])
    #print(lmin==lmin2, lmin, lmin2);



if(args.plot):
    import matplotlib.pyplot as plt 
    
    binding_width = [len(x) for x in gcdrops_regions if x[6] == '1']
    binding_depth = [float(x.score) for x in gcdrops_regions if x[6] == '1']
    nonbinding_width = [len(x) for x in gcdrops_regions if x[6] == '0']
    nonbinding_depth = [float(x.score) for x in gcdrops_regions if x[6] == '0']
    
    plt.plot(binding_depth, binding_width, 'b.', label = 'Bound GC drops')
    plt.plot(nonbinding_depth, nonbinding_width, 'r.', label = 'Unbound GC drops')
    plt.xlabel('minimum GC content')
    plt.ylabel('length of GC drop')
    plt.gca().invert_xaxis()
    print(len(binding_width))
    plt.legend()
    plt.show()
    
    
    
    
    #selected = [x for x in regions if float(x.score) > 500]
    #selected = [(x[0]-500, x[0]+200) for x in enumerate(minima) if x[1]<-7.0]
    #selected = [(x[0]-300, x[0]+300) for x in enumerate(gcmeans) if x[1]<0.142]
    #selected = [(x, x+600) for x in range(0, len(coverage)-300, 300) if  0.22<min(gcmeans[x:x+600])< 0.24]
    #print(selected)

    #flank = 150;
    
    #for region in selected:
        #lcov = coverage[region[0]: region[1]] 
        #lgc = gc_content[region[0]: region[1]] 
        ##lcov = coverage[region.start-flank:region.end+flank]
        ##lgc = gc_content[region.start-flank:region.end+flank]

        #fig, ax1 = plt.subplots()

        #ax1.plot(lcov, 'b-')
        #ax1.set_xlabel("position (nt)")
        #ax1.set_ylabel('coverage', color='b')
        #ax1.tick_params('y', colors='b')
        
        #ax2 = ax1.twinx()
        #ax2.plot(lgc, 'r-')
        #ax2.set_ylabel("gc", color='r')
        #ax2.tick_params('y', colors='r')
        
        #ax1.spines['right'].set_visible(False)
        #ax1.spines['top'].set_visible(False) 

        #fig.tight_layout()
        #plt.show()



