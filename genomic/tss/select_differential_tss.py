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
#parser.add_argument('--outdir', required = True, nargs = '?', type = str, help = "Path to the output plot folder");
#parser.add_argument('--format', default = 'png', nargs = '?', type = str, help = "Plot format");
args = parser.parse_args()


def check(labels, var_names, expression, variants , minexpr, mindiff):
    expr_passed = [];
    for c1, c2 in combinations(range(len(expression)), 2):
        if(min(expression[c1])>=minexpr and min(expression[c2])>=minexpr):
            expr_passed.append((c1,c2))
    result = [];
    for c1, c2 in expr_passed:
        for v in range(len(variants)):
            var1, var2 = variants[v][c1], variants[v][c2]
            diff = max(( min(var1) - max(var2), min(var2) - max(var1) ))
            if(diff > mindiff):
                result.append(( labels[c1], labels[c2], var_names[v], expression[c1], expression[c2], var1, var2 ))
                print(result[-1])
    if(result):
        print()
                
                

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
        
        

for path in [x for x in get_only_files(args.path) if x.endswith("tsv")]:
    labels, var_names, expression, variants = parse_file(path);
    check(labels, var_names, expression, variants , args.minexpr, args.mindiff)
