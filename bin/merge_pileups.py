#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Merges genomic coverage for bed-pileup files'''

import argparse
import os
import sys
from collections import defaultdict, namedtuple
import math

import numpy as np;
import pandas as pd;


parser = argparse.ArgumentParser(description='Merges genomic coverage for bed-pileup files');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the genomic coverage, bed3 format");
#parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Output destination for the statistics plots");
args = parser.parse_args();

total = []

for path in args.path:
    if(len(total)):
        total = total + pd.read_csv(path, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values
    else:
        total = pd.read_csv(path, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values
        


table = pd.read_csv(args.path[0], sep="\t", names = ["chr", "position", "coverage"])
for row, t in zip(table.itertuples(index=False), total):
    print("%s\t%d\t%d" % (row[0], row[1], t))