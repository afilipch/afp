#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Splits transcript according their origin to phage and renormalize TPMs separately for phage and non-phage transcripts'''

import argparse
import os
import sys
from collections import defaultdict
import copy


import numpy as np;
from pybedtools import BedTool, Interval
from itertools import combinations, permutations
import matplotlib.pyplot as plt;




parser = argparse.ArgumentParser(description='Splits transcript according their origin to phage and renormalize TPMs separately for phage and non-phage transcripts');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the expression transcripts, gff format");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args()

def renormalize(transcripts, path, labels):
    expr_list = [];
    for transcript in transcripts:
        expr_list.append([])
        shape = [0];
        for s in transcript.attrs['expression'].split(":"):
            expr_list[-1].extend([float(x) for x in s.split(",")])
            shape.append(len(expr_list[-1]))
    

    expr_list = np.array(expr_list)
    total_expr = np.sum(expr_list, axis=0)
    norma = 1000000/total_expr
    total_expr = [np.mean(total_expr[x[0]:x[1]]) for x in zip(shape, shape[1:])]
    
    with open(path, 'w') as f:
        f.write("#labels=%s\n" % ",".join(labels))
        total_expr = ["%1.3f" % x for x in total_expr]
        f.write("#total_expr=%s\n" % ",".join(total_expr))
        for transcript in transcripts:
            #print(transcript.attrs['expression'].split(":"))
            a = []
            pos = 0;
            for s in transcript.attrs['expression'].split(":"):
                new_expr = [float(x[1])*norma[x[0]+pos] for x in enumerate(s.split(",")) ]
                pos += len(new_expr)
                a.append(",".join(["%1.3f" % x for x in new_expr]))
            #print(a)
            #print()
            transcript.attrs['expression'] = ":".join(a)
            f.write(str(transcript))
        
            
        
    
    


with open(args.path) as f:
    labels = next(f).strip().split("=")[1].split(",")

transcripts = BedTool(args.path)
phage_lists = [x for x in transcripts if int(x.attrs['phage'])], [x for x in transcripts if not int(x.attrs['phage'])]
#print([len(x) for x in phage_lists])
ph_names = ["phage.gff", "nonphage.gff"]

for transcripts, ph_name in zip(phage_lists, ph_names):
    path = os.path.join(args.outdir, ph_name)
    renormalize(transcripts, path, labels)
    

