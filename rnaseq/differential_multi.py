#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Compares multiple rna-seq samples (with/wo replicates) and assigns gene differential expression for all the pairwise combinations of the samples'''

import argparse
import os
import sys
from collections import defaultdict
import copy


import numpy as np;
from scipy.stats import variation
from pybedtools import BedTool, Interval
from itertools import combinations, permutations
#import matplotlib.pyplot as plt;




parser = argparse.ArgumentParser(description='Compares multiple rna-seq samples (with/wo replicates) and assigns gene differential expression for all the pairwise combinations of the samples');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the expression table, tsv format");
parser.add_argument('--minexpr', nargs = '?', default=20, type = float, help = "Minimum mean (among replicates) coverage for a gene to be considered as expressed greatly than the other");
parser.add_argument('--minfold', nargs = '?', default=2.0, type = float, help = "Minimum fold difference between two genes to be considered as differential");
parser.add_argument('--maxfold', nargs = '?', default=20, type = float, help = "An arbitrary number which is set to avoid the division by zero expression");
parser.add_argument('--foldmethod', nargs = '?', default='minmax', choices = ['mean', 'minmax'], type = str, help = "Fold calculation method, either mean1/mean2 [mean] or min1/max2 (if 1>2) [minmax]");

parser.add_argument('--scoring', nargs = '?', default='conservative', choices = ['conservative', 'explorative'], type = str, help = "Scoring type with more focus on robustness or interesting candidates");
args = parser.parse_args()


def scoring_conservative(coeff, e):
    return coeff*e

def scoring_explorative(coeff, e):
    return coeff*(e**(0.5))#(np.log2(e+2))


if(args.scoring == 'conservative'):
    scoring = scoring_conservative;
if(args.scoring == 'explorative'):
    scoring = scoring_explorative;
    

if(args.foldmethod == 'mean'):
    def get_fold_change(expr1, expr2):
        return min(expr1)/max(expr2)
elif(args.foldmethod == 'minmax'):
    def get_fold_change(expr1, expr2):
        return np.mean(expr1)/np.mean(expr2)


def gene_total_score(mylist):
    return max([x[1] for x in mylist]), np.argmax([x[1] for x in mylist])
    
    
    
    

def assign_score(expr1, expr2, fold, change):
    e = max([min(expr1), min(expr2)])
    if(not change):
        coeff = 0.5*(fold - 1)
    else:
        coeff = fold;
        
    return scoring(coeff, e)







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



def compare(expr1, expr2, minfold, minexpr, maxfold):
    
    if(any(expr1)):
        if(any(expr2)):
            f1 = get_fold_change(expr1, expr2);
            f2 = get_fold_change(expr2, expr1)
            change, fold = check_change(expr1, expr2, f1, f2, minfold, minexpr)
        else:
            fold, change = maxfold, simple_check_change(expr1, minexpr)
    else:
        if(any(expr2)):
            fold, change = 1/maxfold, 2*simple_check_change(expr2, minexpr)
        else:
            fold, change = 1, 0
            
    if(change == 2):
        fold = 1/float(fold)
        
    return fold, change
        

def get_expression_line(transcript):
    e_line = transcript.attrs['expression'].split(":")
    e_arr = [[float(y) for y in x.split(",")] for x in e_line]
    mean_line = ["%1.2f" % np.mean(x) for x in e_arr]
    var_line = ["%1.1f" % (variation(x)*100) if any(x) else '0' for x in e_arr]
    return e_line + mean_line + var_line
        
        
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
        fold, change = compare(expr1, expr2, args.minfold, args.minexpr, args.maxfold)
        score = assign_score(expr1, expr2, fold, change)
        fold_list[c].append((fold, score, change))
        
 
header = ['gene', 'pair', 'score', 'change'] + labels + ["%s (mean)" % x for x in labels] + ["%s (variation)" % x for x in labels] + ["|".join(x) for x in labels_combinations]
print("\t".join(header))
for c, transcript in enumerate(transcripts):
    fold_list_local = [x[c] for x in fold_list]
    ##print(len(fold_list_local), len(labels_combinations))
    score, pos = gene_total_score(fold_list_local)
    label_pair = labels_combinations[pos]
    change = fold_list_local[pos][2]
    a = [transcript.attrs['ID'], "|".join(label_pair), "%d" % score, str(change)] + get_expression_line(transcript) + ["%1.2f" % x[0] for x in fold_list_local]
    print("\t".join(a))
    

