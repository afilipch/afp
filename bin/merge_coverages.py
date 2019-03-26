#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Merges provided coverage tracks. That is, sums up coverage values for each position'''

import argparse
import os
import sys
import numpy as np;
import pandas as pd;


parser = argparse.ArgumentParser(description='Merges provided coverage tracks. That is, sums up coverage values for each position');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the coverage files");
args = parser.parse_args();

###Read Coverage
data = pd.read_csv(args.path[0], sep="\t" , names = ["chr", "position", "coverage"])
coverage = data.coverage.values

for path in args.path[1:]:
    coverage += pd.read_csv(path, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values

data.coverage = coverage;
for l in data.itertuples():
    print("%s\t%s\t%d" % (l.chr, l.position, l.coverage));    

