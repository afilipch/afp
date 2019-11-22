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

data = pd.read_csv(args.coverage, sep="\t" , names = ["chr", "position", "coverage"])
coverage = copy.copy(data.coverage.values)


peaks = [x for x in BedTool(args.path) if float(x.score) > args.zscore]
peak_coverage = []

for peak in peaks:
    peak_coverage.extend(coverage[peak.start:peak.end]) 
 
sys.stderr.write("###normalize_coverage\n")
sys.stderr.write("%s\nFile %s is processed\n\nSignal fraction of the coverage:%1.3f\t\nNoise fraction of the coverage:%1.3f\t\n\n" % ("_"*120, args.path, sum(peak_coverage)/sum(coverage), 1 - sum(peak_coverage)/sum(coverage)))

normfactors = {}
normfactors['all'] = np.mean(coverage);
normfactors['signal'] = np.mean(peak_coverage);
normfactors['noise']= (sum(coverage) - sum(peak_coverage))/(len(coverage) - len(peak_coverage))


data.coverage = coverage/normfactors[args.mode];
all_coverage = coverage/normfactors['all']
signal_coverage = coverage/normfactors['signal']
noise_coverage = coverage/normfactors['noise']

for l, ac, sc, nc, oc in zip(data.itertuples(), all_coverage, signal_coverage, noise_coverage, coverage):
    print("%s\t%s\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%1.5f" % (l.chr, l.position, l.coverage, ac, sc, nc, oc));

