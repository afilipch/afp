#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Normalizes coverage track for chip/chap-seq experiments in order to compare different samples'''

import argparse
import os
import sys
import copy
from collections import defaultdict


import numpy as np;
import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Normalizes coverage track for chip/chap-seq experiments in order to compare different samples');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to all the detected peaks");
parser.add_argument('--coverage', required = '?', default=True, type = str, help = "Path to the coverage files, must have the order consistent with input \'peak\' files")
parser.add_argument('--zscore', nargs = '?', default=2.0, type = float, help = "Mininimum z-score required for a peak to be considered as signal");

parser.add_argument('--mode', nargs = '?', choices=['all', 'noise', 'signal'], default='all', type = str, help = "Mode of the normalization. --all: normalize to the whole coverage, --signal: noramlize to the coverage of the filtered (with z-score) peaks, --noise: normalize to the coverage of the non-signal regions");
args = parser.parse_args();

table = pd.read_csv(args.coverage, sep="\t" , names = ["chr", "postion", "coverage"])
coverage = table.coverage.values


peaks = [x for x in BedTool(args.path) if float(x.score) > args.zscore]
length = len(coverage)
temp = copy.copy(coverage)

for peak in peaks:
    temp[peak.start:peak.end] = [0]*len(peak)
    length = length - len(peak);
    
sys.stderr.write("%s\nFile %s is processed\n\nSignal fraction of the coverage:%1.3f\t\nNoise fraction of the coverage:%1.3f\t\n\n" % ("_"*120, args.path, 1-sum(temp)/sum(coverage), sum(temp)/sum(coverage)))


if(args.mode == 'all'):
    normfactor = np.mean(coverage);

elif(args.mode == 'noise'):
    normfactor = sum(temp)/length
    
elif(args.mode == 'signal'):
    temp = []
    for peak in peaks:
        temp.extend(coverage[peak.start:peak.end])
    normfactor = np.mean(temp);


coverage_normed = coverage/normfactor
for cov, row in zip(coverage_normed, table.itertuples()):
    print("\t".join((row[1], str(row[2]), "%1.5f" % cov)))

