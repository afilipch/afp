#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts tss from coverage to bed format'''
import argparse
import os
import sys
import numpy as np;
import matplotlib.pyplot as plt;
from collections import defaultdict

from pybedtools import BedTool
from Bio import SeqIO

from afbio.sequencetools import coverage2dict

parser = argparse.ArgumentParser(description='Converts tss from coverage to bed format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the TSS in coverage format");
parser.add_argument('--threshold', nargs = '?', default=10, type = float, help = "Threshold (in medians) for the TSS intensity");
parser.add_argument('--strand', nargs = '?', required=True, choices = ['+', '-'], type = str, help = "Strand of the provided coverage");
args = parser.parse_args();



def stat_cov(covdict):
    allcov  = [];
    for v in covdict.values():
        allcov.extend(v)
    ###Output basic coverage statistics
    sys.stderr.write("###Explore TSS\n")
    signal = [x for x in allcov if x]
    median = np.median(signal)
    sys.stderr.write( "Number of TSS:\t%d\nMedian TSS intensity:\t%1.2f\nMean TSS intensity:\t%1.2f\nMax TSS intensity:\t%1.2f\n\n" % (len(signal), median, np.mean(signal), max(allcov)) )
    return median

covdict = coverage2dict(args.path)
median = stat_cov(covdict)
threshold = args.threshold*median
passed, out = 0, 0

for chrom, coverage in covdict.items():
    for pos, intensity in enumerate(coverage):
        if(intensity):
            if(intensity >= threshold):
                print("%s\t%d\t%d\t.\t%d\t%s" % (chrom, pos, pos+1, intensity, args.strand));
                passed += 1;
            else:
                out += 1;
                
sys.stderr.write("Passed TSS:\t%d\nFiltered out TSS:\t%d\n\n" % (passed, out) )


