#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates table with all the information regarding HrrA project'''

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




parser = argparse.ArgumentParser(description='Creates table with all the information regarding HrrA project');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated consensus regions of the time-series binding peaks, '~/afp/postanalyses/reannotate_according_tss.py' output");
parser.add_argument('--diff', nargs = 3, required=True, type = str, help = "Path to the file with genes\' differential expression");
parser.add_argument('--chapnames', nargs = '+', required=True, type = str, help = "Names for CHAP experiment");
parser.add_argument('--diffnames', nargs = 3, required=True, type = str, help = "Names for diff RNA-Seq experiment");
#parser.add_argument('--distance', nargs = '?', default=400, type = int, help = "Maximum allowed distance for the peak to the closest peak");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Path to the plot directory");
args = parser.parse_args();


header = ["Chrom", "Start", "Stop", "Gene ID", "Gene symbol", "Distance ATG", "Distance to TSS"] + ["ChAP T=%s" % x for x in args.chapnames] + ["mRNA wt T=%s" % x for x in args.diffnames] + ["mRNA DhrrA T=%s" % x for x in args.diffnames] + ["Log2 DhrrA/WT T=%s" % x for x in args.diffnames] + ["Predicted Function", "Annotation"]

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
            
            gene2diff[name].append((wt, ko, log_ko_wt))
            
            
regions = BedTool(args.path);
#sys.stderr.write(str(len(regions)) + "\n\n")

print("\t".join(header));
ll = [];
for region in regions:
    curdiff = gene2diff[region.attrs['gene']]
    if(len(curdiff) == 3):
        a = [region.chrom, str(region.start), str(region.stop), region.attrs['gene'],  region.attrs['genesymbol'], region.attrs['atg'], region.attrs['tss']] + ["%1.3f" % float(x) if x != "None" else "None" for x in region.attrs['maxcov'].split(",")] + ["%1.1f" % x[0] for x in curdiff] + ["%1.1f" % x[1] for x in curdiff] + ["%1.3f" % x[2] for x in curdiff] + [region.attrs['function'],  region.attrs['annotation']]
        print("\t".join(a))
    #else:
        #print(region.attrs['genesymbol'], curdiff)
    ll.append(len(curdiff))

sys.stderr.write( "%s\n" % Counter(ll))
    
            
            
            
#for k, v in gene2diff.items():
    #print(len(v));
            
            
