#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Filters the detected peaks based on the distribution of their scores'''

import argparse
import os
import sys
import scipy
import numpy as np;
import pandas as pd;
import matplotlib.pyplot as plt;
import os
from pybedtools import BedTool
from collections import defaultdict
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Filters the detected peaks based on the distribution of their scores');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the detected peaks with assigned zscore");
parser.add_argument('--coverage', nargs = '?', required=True, type = str, help = "Path to the coverage track");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
parser.add_argument('--convolution', nargs = '?', required=True, type = str, help = "Path to the convolved coverage track");
parser.add_argument('--flen', nargs = '?', default=200, type = int, help = "Length of the peak\'s flanks to be drawn");
parser.add_argument('--plot', nargs = '?', default=False, const=True, type = bool, help = "");
args = parser.parse_args();

coverage = pd.read_csv(args.coverage, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values
convolution = pd.read_csv(args.convolution, sep="\t" , names = ["chr", "position", "convolution"]).convolution.values
peaks = BedTool(args.path);
zscores = [float(x.score) for x in peaks]
peaks = [x for x in peaks if float(x.score)>2]
#print(zscores)

zlimits = [min(zscores)-1, 0, 2, 3, 4, 5, 7, 10, 20, 50]
limitnames = ["_%d" % zlimits[1]] + ["%d_%d" % x for x in zip(zlimits[1:-1], zlimits[2:])] + ["%d_" % zlimits[-1]]


def definerange(score, zlimits):
    for c, zl in enumerate(zlimits[::-1]):
        if(score>zl):
            return(len(zlimits) - c - 1);
    else:
        return 0;
    

peaksdict = defaultdict(list);
for peak in peaks:
    peaksdict[limitnames[definerange(float(peak.score), zlimits)]].append(peak);


def makefig(l_coverage, l_convolution, output):
    fig, ax1 = plt.subplots()

    ax1.plot(l_coverage, 'b-')
    ax1.set_xlabel("postion (nt)")
    ax1.set_ylabel('coverage', color='b')
    ax1.tick_params('y', colors='b')
    
    l_convolution = [max(0,x) for x in l_convolution]
    ax2 = ax1.twinx()
    ax2.plot(l_convolution, 'r-')
    ax2.set_ylabel("convolution", color='r')
    ax2.tick_params('y', colors='r')

    fig.tight_layout() 
    plt.savefig(output, format='svg')   
    plt.close()



#print(peaksdict)

###Plot coverage vs convolution
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
    
for zln, plist in peaksdict.items():
    limitdir = os.path.join(args.outdir, zln)
    if not os.path.exists(limitdir):
        os.makedirs(limitdir)
    
    for c, peak in enumerate(plist):
        start = max(0, peak.start - args.flen)
        end = peak.stop + args.flen
        l_coverage = coverage[start:end]
        l_convolution = convolution[start:end]
        makefig(l_coverage, l_convolution, os.path.join(limitdir, "peak_%d.svg" % (c+1)))
    
#if(args.plot):

    #fig, ax1 = plt.subplots()

    #ax1.plot(coverage, 'b-')
    #ax1.set_xlabel("postion (nt)")
    #ax1.set_ylabel('coverage', color='b')
    #ax1.tick_params('y', colors='b')
    
    #convolution = [max(0,x) for x in convolution]
    #ax2 = ax1.twinx()
    #ax2.plot(convolution, 'r-')
    #ax2.set_ylabel("convolution", color='r')
    #ax2.tick_params('y', colors='r')

    #fig.tight_layout() 
    #plt.show()
    

    
#print(limitnames[definerange(0.7, zlimits)])    

#print(limitnames)
#
