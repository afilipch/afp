#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts consensus binding regions into more human-readable tsv file'''

import argparse
import sys
import os
from collections import defaultdict
from bisect import bisect_right, bisect_left


import pandas as pd;
import numpy as np;
import matplotlib.pyplot as plt;
from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Converts consensus binding regions into more human-readable tsv file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the consensus regions, gff format");
parser.add_argument('--names', nargs = '+', required=True, type = str, help = "Names of the samples. Must be provided in the order consistent with regions");
args = parser.parse_args();

# chr1    un      consensus_region        195051  195131  0       +       .       Name=195090;peakpos=195090,195110,195095,195104,195093,195115;zscores=9.5,16.1,25.7,25.1,17.0,14.7;maxcov=6.721,7.178,9.741,9.855,7.827,7.865;start_gene=hkm;start_gene_distance=52;start_gene_strand=-;start_annotation=putative sensor histidine kinase fragment of two-component system, putative pseudogene;start_function=Post-translational modification;

print("\t".join(('chr', 'start', 'end', 'top_position', 'start_gene', 'start_gene_distance', 'start_gene_strand',  'start_function', 'start_annotation', "\t".join(["zscore_%s" % x for x in args.names]), "\t".join(["intensity_%s" % x for x in args.names]))))
for interval in BedTool(args.path):
    print("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (interval.chrom, interval.start, interval.end, interval.name, interval.attrs['start_gene'], interval.attrs['start_gene_distance'], interval.attrs['start_gene_strand'] ,interval.attrs['start_function'], interval.attrs['start_annotation'], "\t".join(['0' if x == 'None' else x for x in interval.attrs['zscores'].split(",")]),  "\t".join(['0' if x == 'None' else x for x in interval.attrs['maxcov'].split(",")])   ))