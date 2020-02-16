#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Selects the genes which have differential tss expression'''

import argparse
import os
import sys
from collections import defaultdict
import copy


import numpy as np;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;
from itertools import combinations, permutations

from afbio.generators import get_only_files


parser = argparse.ArgumentParser(description='Selects the genes which have differential tss expression');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with tss expression files, tsv format");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts, gff format");
parser.add_argument('--minexpr', nargs = '?', default = 100, type = float,  help = "Minimum allowed expression (TPM)");
parser.add_argument('--mindiff', nargs = '?', default = 0.25, type = float,  help = "Minimum allowed difference in TSS fractions (normalized to 1)");
parser.add_argument('--minfraction', nargs = '?', default = 2, type = float,  help = "Minimum allowed fraction to assign a spurious transcription");
parser.add_argument('--outdir', required = True, nargs = '?', type = str, help = "Path to the output folder");
#parser.add_argument('--format', default = 'png', nargs = '?', type = str, help = "Plot format");
args = parser.parse_args()

def get_score(e1, e2, diff):
    return (min(e1) + min(e2))*diff


#def check_spurious_transcription(var, minfraction):
    #if(min(var) > minfraction)
    
    

def check(labels, var_names, expression, variants , minexpr, mindiff, minfraction):
    expr_passed = [];
    for c1, c2 in combinations(range(len(expression)), 2):
        if(min(expression[c1])>=minexpr and min(expression[c2])>=minexpr):
            expr_passed.append((c1,c2))
            
    normal = [];
    spurious = [];
    for c1, c2 in expr_passed:
        for v in range(len(variants)):
            var1, var2 = variants[v][c1], variants[v][c2]
            if(min(var1) > minfraction):
                score = min(var1)*min(expression[c1])
                spurious.append(( labels[c1], var_names[v], expression[c1], var1, score))
            
            diff = max(( min(var1) - max(var2), min(var2) - max(var1) ))
            if(diff > mindiff):
                score = get_score(expression[c1], expression[c2], diff)
                normal.append(( labels[c1], labels[c2], var_names[v], expression[c1], expression[c2], var1, var2, score))
    if(spurious):
        #print(max(spurious, key = lambda x: x[-1]))
        return [], max(spurious, key = lambda x: x[-1])
    elif(normal):
        return max(normal, key = lambda x: x[-1]), []
    else:
        return [], []
                
                

def parse_file(path):
    with open(path) as f:
        labels = next(f).strip().split("\t")[1:]
        expression = [[float(y) for y in x.split(",")] for x in next(f).strip().split("\t")[1:]]
        variants, var_names = [], [];
        for l in f:
            a = l.strip().split("\t")
            variants.append([ [float(y) for y in x.split(",")] for x in a[1:]])
            var_names.append(a[0])
        return labels, var_names, expression, variants
        
        
normal_list, spurious_list = [], []
for path in [x for x in get_only_files(args.path) if x.endswith("tsv")]:
    name = os.path.basename(path).split(".")[0]
    labels, var_names, expression, variants = parse_file(path);
    res =  check(labels, var_names, expression, variants , args.minexpr, args.mindiff, args.minfraction)
    if(res[0]):
        normal_list.append((name, res[0]))
    if(res[1]):
        spurious_list.append((name, res[1]))
        
        
normal_list.sort(reverse=True, key = lambda x: x[1][-1])
spurious_list.sort(reverse=True, key = lambda x: x[1][-1])        
with open(os.path.join(args.outdir, "normal.tsv"), 'w') as f:
    header = ["gene", "time1", "time2", "tss", "score", "expression1", "expression2", "fraction1", "fraction2"]
    f.write("%s\n" % "\t".join(header) )
    for name, l in normal_list:
        a = [",".join([str(y) for y in x]) for x in l[3:7]]
        f.write("%s\n" % "\t".join([name] + list(l[:3]) + [str(int((l[-1])))] + a) )
        
        
with open(os.path.join(args.outdir, "spurious.tsv"), 'w') as f:
    header = ["gene", "time", "tss", "score", "expression", "fraction"]
    f.write("%s\n" % "\t".join(header) )
    for name, l in spurious_list:
        a = [",".join([str(y) for y in x]) for x in l[2:4]]
        f.write("%s\n" % "\t".join([name] + list(l[:2]) + [str(int((l[-1])))] + a) )
        
        
        
        
        
        
