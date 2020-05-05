#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Extracts genomic coverage for the provided regions'''

import argparse
import os
import sys
from collections import defaultdict

import numpy as np;
import matplotlib.pyplot as plt;

from pybedtools import BedTool

from afbio.sequencetools import coverage2dict


parser = argparse.ArgumentParser(description='Extracts genomic coverage for the provided regions');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "Path to the coverage files (plus and minus), bed format");
parser.add_argument('--regions', nargs = '?', required=True, type = str, help = "Path to the regions of interest, gff/bed format");
args = parser.parse_args();

#ORDER = ['CDS', 'rRNA', 'tRNA', 'pseudogene']

covdict = {'+': coverage2dict(args.path[0]), '-': coverage2dict(args.path[1])}

for interval in BedTool(args.regions):
    local_cov =  covdict[interval.strand][interval.chrom][interval.start:interval.stop]
    for pos, el in enumerate(local_cov):
        print(interval.name, pos, el)
        
    
