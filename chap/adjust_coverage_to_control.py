#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Adjusts genomic coverage based on the provided control (dna expresssion) coverage'''

import argparse
import os
import sys
import numpy as np;
import pandas as pd;


parser = argparse.ArgumentParser(description='Adjusts genomic coverage based on the provided control (dna expresssion) coverage');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the coverage file, bed3 format");
parser.add_argument('--control', nargs = '?', required=True, type = str, help = "Path to the control file, bed3 format");
args = parser.parse_args();

###Read Coverage
data = pd.read_csv(args.path, sep="\t" , names = ["chr", "position", "coverage"])
coverage = data.coverage.values


control_coverage = pd.read_csv(args.control, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values
control_coverage += 1;
control_coverage = control_coverage/np.mean(control_coverage)


data.coverage = coverage/control_coverage;
for l in data.itertuples():
    print("%s\t%s\t%1.1f" % (l.chr, l.position, l.coverage)); 
