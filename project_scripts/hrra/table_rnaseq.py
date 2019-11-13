#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates rna-seq table for all the genes'''

import argparse
import os
import sys
#import copy
from collections import defaultdict, Counter
from itertools import product;

import numpy as np;
from scipy.stats import pearsonr
#import pandas as pd;
from pybedtools import BedTool, Interval




parser = argparse.ArgumentParser(description='Creates rna-seq table for all the genes');
parser.add_argument('--diff', nargs = 3, required=True, type = str, help = "Path to the file with genes\' differential expression");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the geneid/genesymbol connection table");
parser.add_argument('--diffnames', nargs = 3, required=True, type = str, help = "Names for diff RNA-Seq experiment");

args = parser.parse_args();


header = ["Gene ID", "Gene symbol"] + ["mRNA wt T=%s" % x for x in args.diffnames] + ["mRNA DhrrA T=%s" % x for x in args.diffnames] + ["Log2 DhrrA/WT T=%s" % x for x in args.diffnames] + ["Predicted Function", "Annotation"]
print("\t".join(header))


ann_dict = {}
cg_dict = {};
with open(args.annotation) as f:
    next(f)
    for l in f:
        a = l.strip().split(";")
        if(len(a) > 9 and a[1] and a[2]):
            ann_dict[a[1]] = a[2], a[8].replace('\"', ''), a[9].replace('\"', '')


gene2diff = defaultdict(list)
for path in args.diff:
    with open(path) as f:
        next(f)
        for l in f:
            a = l.strip().split("\t");
            name = a[0].split("-")[1]
            wt = np.mean([float(x) for x in a[1].split(";")])
            ko = np.mean([float(x) for x in a[2].split(";")])
            log_ko_wt = float(a[3]);
            significant_change = float(a[5]);
            
            gene2diff[name].append((wt, ko, log_ko_wt, significant_change))
            
            
for geneid, curdiff in sorted(gene2diff.items(), key=lambda x: x[0]):
    if( [x for x in curdiff if x[3] and abs(x[2])>=1] ):
        genesymbol, annotation, function = ann_dict.get(geneid, ('unknown', 'unknown', 'unknown'))
        a = [geneid, genesymbol] + ["%1.1f" % x[0] for x in curdiff] + ["%1.1f" % x[1] for x in curdiff] + ["%1.3f" % x[2] for x in curdiff] + [function,  annotation] 
        print("\t".join(a))
        
