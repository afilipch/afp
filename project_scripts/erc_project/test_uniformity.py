#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Checks whether prophages are distributed not uniformly along host genome'''

import argparse
import os
import sys

import numpy as np;

from pybedtools import BedTool, Interval
from scipy.stats import binom_test




parser = argparse.ArgumentParser(description='Checks whether prophages are distributed not uniformly along host genome');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the prophage table, gff format");
parser.add_argument('--maxlength', nargs = '?', default=10e10, type = int, help = "Maximal length allowed for prophages");
parser.add_argument('--btype', nargs = '?', default='Bacteria', type = str, help = "Type of bacteria to be analyzed");
parser.add_argument('--maxfraction', nargs = '?', default=0.05, type = float, help = "Maximal allowed prophage genome length (as a fraction of host genome length) to be reported as circle");
args = parser.parse_args();

print("%s\n" % args.btype)


### Reading the input

prophages = BedTool(args.path)

# filter length
prophages_filtered = [x for x in prophages if len(x) <= args.maxlength]
sys.stderr.write("%d prophages passed the length filtering (length <= %d) out of %d\n\n" % (len(prophages_filtered), args.maxlength, len(prophages)))
prophages = prophages_filtered;
prophages_filtered = [x for x in prophages if args.btype in x.attrs['host_family']]
sys.stderr.write("%d %s prophages are selected from total %d prophages\n\n" % (len(prophages_filtered), args.btype, len(prophages)))
prophages = prophages_filtered;


ptypes = ['viable', 'remnant']
for integrase, ptype in zip( ('True', 'False'), ptypes):
    
    local_prophages = [x for x in prophages if x.attrs['integrase'] == integrase]
    host_rel_positions = [ (int(x.attrs['host_end']) + int(x.attrs['host_start']))/2 for x in local_prophages]
    bleft, bright = 250, 750
    down = len([x for x in host_rel_positions if (x>=bleft and x<bright)])
    total = len(host_rel_positions)
    pval = binom_test(down, n=total, p=0.5)

    print("%s prophages: %d\nNumber of prophages in the half of bacterial chromosome opposite to replication origin %d\nBinomial test p-value %s\n" % (ptype, total, down, '{:.2e}'.format(pval)))

print("\n%s" % ("_"*100) )
