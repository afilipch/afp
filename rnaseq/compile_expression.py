#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Compiles gene expression values coming form different sources'''

import argparse
import os
import sys
from collections import defaultdict

import numpy as np;
import matplotlib.pyplot as plt;

from pybedtools import BedTool
from afbio.generators import get_only_files


parser = argparse.ArgumentParser(description='Compiles gene expression values coming form different sources');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the transcripts folder");
#parser.add_argument('--first', nargs = '+', required=True, type = str, help = "Path to gene expression data (replicates) for the first sample");
#parser.add_argument('--second', nargs = '+', required=True, type = str, help = "Path to gene expression data (replicates) for the second sample");
#parser.add_argument('--labels', nargs = 2, required=True, type = str, help = "Names of the samples");
parser.add_argument('--order', nargs = '+', type = str, help = "If set, labels are output in this order");
args = parser.parse_args();


def get_tpms(path):
    return dict([(x.name, float(x.attrs['tpm'])) for x in BedTool(path)])

def line2score(l):
    return sum([sum(x) for x in l[1:]])

    
    

label2file = defaultdict(list)
for f in sorted(get_only_files(args.path)):
    label = "_".join(os.path.basename(f).split("_")[:-1])
    label2file[label].append(f)
    
if(args.order):
    label2file = [(x, label2file[x]) for x in args.order]
else:
    label2file = list(sorted(label2file.items(), key = lambda x: x[0]))
  
expression = [];
labels = [];
genes = set()
for label, local_files in label2file:
    local_expr = [];
    for lf in local_files:
        #print(lf)
        gd = get_tpms(lf)
        #print(gd)
        local_expr.append(gd)
        genes.update(gd.keys());
    expression.append(local_expr)
    labels.append(label)
    
    
    
  
#print(genes)
table_list = [];
for gene in genes:
    line_list = [gene];
    for local_expr in expression:
        line_list.append([x[gene] for x in local_expr])
    table_list.append(line_list)

print("\t".join(['gene'] + labels))    
table_list.sort(key = lambda x: line2score(x), reverse = True)
for l in table_list[:10]:
    s_list = []
    for a in l[1:]:
        s_list.append([str(x) for x in a])
        
    s = "\t".join([l[0]] + [",".join(x) for x in s_list])
    print(s)
    


sys.exit();
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




