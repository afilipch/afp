#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Compares multiple rna-seq samples (with/wo replicates) and assigns gene differential expression for all the pairwise combinations of the samples'''

import argparse
import os
import sys
from collections import defaultdict
import copy


import numpy as np;
from pybedtools import BedTool, Interval
from itertools import combinations, permutations
import matplotlib.pyplot as plt;
import yaml;



parser = argparse.ArgumentParser(description='Compares multiple rna-seq samples (with/wo replicates) and assigns gene differential expression for all the pairwise combinations of the samples');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the expression table, tsv format");
parser.add_argument('--minexpr', nargs = '?', default=20, type = float, help = "Minimum mean (among replicates) coverage for a peak to be considered as expressed greatly than the other");
parser.add_argument('--minfold', nargs = '?', default=2.0, type = float, help = "Minimum fold difference between two peaks to be considered as differential");
#parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args()



def gene_total_score(mylist):
    return max([x[1] for x in mylist]), np.argmax([x[1] for x in mylist])
    
    
    
    

def assign_score(expr1, expr2, fold, change, maxfold=20):
    e = max([min(expr1), min(expr2)])
    
    if (fold == 'NaN'):
        myfold = maxfold;
    else:
        myfold = min(maxfold, fold)
        
    
    if(not change):
        coeff = 0.5*(myfold - 1)
    else:
        coeff = myfold;
        
    return coeff*e




def get_fold_change(expr1, expr2):
    return min(expr1)/max(expr2)



def simple_check_change(expr, minexpr):
    return int(all([x>minexpr for x in expr]))


def check_change(expr1, expr2, f1, f2, minfold, minexpr):
    m1 = max(expr1)
    m2 = max(expr2)
    change = 0;
    
    if(f1>=f2):
        fold = f1;
        if(f1>minfold and all([x>minexpr for x in expr1])):
            change = 1;
    else:
        fold = f2;
        if(f2>minfold and all([x>minexpr for x in expr2])):
            change = 2;        
        
    return change, fold



def compare(expr1, expr2, minfold, minexpr):
    
    if(any(expr1)):
        if(any(expr2)):
            f1 = get_fold_change(expr1, expr2);
            f2 = get_fold_change(expr2, expr1)
            change, fold = check_change(expr1, expr2, f1, f2, minfold, minexpr)
        else:
            fold, change = 'NaN', simple_check_change(expr1, minexpr)
    else:
        if(any(expr2)):
            fold, change = 'NaN', 2*simple_check_change(expr2, minexpr)
        else:
            fold, change = 'NaN', 0
    return fold, change
        
    
        
        
with open(args.path) as f:
    labels = next(f).strip().split("=")[1].split(",")



transcripts = BedTool(args.path)
expr_list = [[] for i in range(len(labels))];

         
for transcript in transcripts:
    for c, s in enumerate(transcript.attrs['expression'].split(":")):
        expr_list[c].append([float(x) for x in s.split(",")])

#print(len(expr_list[0]));
#sys.exit()

labels_combinations = list(combinations(labels, 2))
fold_list = [[] for x in labels_combinations]
for c, (elist1, elist2) in enumerate(combinations(expr_list, 2)):
    for expr1, expr2 in zip(elist1, elist2):
        fold, change = compare(expr1, expr2, args.minfold, args.minexpr)
        score = assign_score(expr1, expr2, fold, change)
        fold_list[c].append((fold, score, change))
        
 
header = ['gene', 'pair', 'score', 'change'] + labels + ["|".join(x) for x in labels_combinations]
print("\t".join(header))
for c, transcript in enumerate(transcripts):
    fold_list_local = [x[c] for x in fold_list]
    ##print(len(fold_list_local), len(labels_combinations))
    score, pos = gene_total_score(fold_list_local)
    label_pair = labels_combinations[pos]
    change = fold_list_local[pos][2]
    a = [transcript.attrs['ID'], "|".join(label_pair), "%d" % score, str(change)] + transcript.attrs['expression'].split(":") + ["%1.2f" % x[0] if x[0] != 'NaN' else 'NaN' for x in fold_list_local]
    print("\t".join(a))
    

