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
import yaml


parser = argparse.ArgumentParser(description='Compiles gene expression values coming form different sources');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the transcripts folder");
#parser.add_argument('--outdir', nargs = '?', default='.', type = str, help = "Path to the output directory");
parser.add_argument('--annotation', nargs = '?', type = str, help = "Path to the gene annotation");
#parser.add_argument('--first', nargs = '+', required=True, type = str, help = "Path to gene expression data (replicates) for the first sample");
#parser.add_argument('--second', nargs = '+', required=True, type = str, help = "Path to gene expression data (replicates) for the second sample");
#parser.add_argument('--labels', nargs = 2, required=True, type = str, help = "Names of the samples");
parser.add_argument('--order', nargs = '+', type = str, help = "If set, labels are output in this order");
args = parser.parse_args();


def get_tpms(path):
    return dict([(x.name, float(x.attrs['tpm'])) for x in BedTool(path)])

def line2score(l):
    return sum([sum(x) for x in l[1:]])

    
gene2annotation = dict([ (x.attrs['ID'], x) for x in BedTool(args.annotation)])


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
    
    
   
   
print("#labels=%s" % ",".join(labels))
table_list.sort(key = lambda x: line2score(x), reverse = True)
for l in table_list:
    transcript = gene2annotation[l[0][5:]]
    s_list = []
    for a in l[1:]:
        s_list.append([str(x) for x in a])
    transcript.attrs['expression'] = ":".join([",".join(x) for x in s_list])
    sys.stdout.write(str(transcript))
    

#with open(os.path.join(args.outdir, 'table.raw.tsv'), 'w') as f:
    #f.write("%s\n" % "\t".join(['gene'] + labels))    
    #table_list.sort(key = lambda x: line2score(x), reverse = True)
    #for l in table_list:
        #s_list = []
        ##print(l)
        #for a in l[1:]:
            #s_list.append([str(x) for x in a])
            
        #s = "\t".join([l[0]] + [",".join(x) for x in s_list])
        #l[0] = list(gene2annotation[l[0][5:]])
        #f.write("%s\n" % s)
    


#with open(os.path.join(args.outdir, 'table.raw.yaml'), 'w') as f:
    #to_yaml = labels, table_list
    #yaml.dump(to_yaml, f)






