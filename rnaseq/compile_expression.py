#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Compiles gene expression values coming form different sources'''

import argparse
import os
import sys
from collections import defaultdict

import numpy as np;
import matplotlib.pyplot as plt;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Compiles gene expression values coming form different sources');
parser.add_argument('--first', nargs = '+', required=True, type = str, help = "Path to gene expression data (replicates) for the first sample");
parser.add_argument('--second', nargs = '+', required=True, type = str, help = "Path to gene expression data (replicates) for the second sample");
parser.add_argument('--labels', nargs = 2, required=True, type = str, help = "Names of the samples");
args = parser.parse_args();


#def get_sampleid(multipath):
    #return os.path.basename(multipath[0]).split(".")[0]

def get_expr(mp1, mp2):
    genes2expr = defaultdict(lambda: defaultdict(list));
    for path in mp1:
        for interval in BedTool(path):
            genes2expr[interval.name][0].append(float(interval.attrs['tpm']))
    for path in mp2:
        for interval in BedTool(path):
            genes2expr[interval.name][1].append(float(interval.attrs['tpm']))
    return genes2expr;
            
genes2expr = [ (x[0], x[1][0], x[1][1]) for x in get_expr(args.first, args.second).items()]
genes2expr.sort(key= lambda x: sum(x[1]), reverse=True)

print("\t".join(["gene_id"] + args.labels));
for el in genes2expr:
    print("\t".join((el[0], ";".join([str(x) for x in el[1]]), ";".join([str(x) for x in el[2]]))))




