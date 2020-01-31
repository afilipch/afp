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
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args()


#with open(args.path) as f:
    #labels, gene_expression_list = yaml.full_load(f)
    
#print(labels)
#print(gene_expression_list[0])

#sys.exit()


def gene_total_score(mylist):
    return max([x[1] for x in mylist]), np.argmax([x[1] for x in mylist])
    
    
    
    

def assign_score(expr1, expr2, fold, change):
    e = max([np.mean(expr1), np.mean(expr2)])
    
    if (fold == 'NaN'):
        myfold = 10;
    else:
        myfold = max(fold, 1/fold)
    
    if(not change):
        coeff = 0.5*(myfold - 1)
    else:
        coeff = 1;
        
    return coeff*myfold*e



def simple_check_change(expr, minexpr):
    return int(all([x>minexpr for x in expr]))


def check_change(expr1, expr2, fold, minfold, minexpr):
    m1 = max(expr1)
    m2 = max(expr2)
    change = 0;
    if(fold>=1 and fold>minfold and all([x>minexpr for x in expr1]) and all([x>m2 for x in expr1])):
        change = 1;
    elif( 1.0/fold>minfold and all([x>minexpr for x in expr2]) and all([x>m1 for x in expr2]) ):
        change = 2;
        
    return change;



def compare(expr1, expr2, minfold, minexpr):
    s1 = np.mean(expr1)
    s2 = np.mean(expr2)
    if(s1):
        if(s2):
            fold = s1/s2
            change = check_change(expr1, expr2, fold, minfold, minexpr)
        else:
            fold, change = 'NaN', simple_check_change(expr1, minexpr)
    else:
        if(s2):
            fold, change = 'NaN', 2*simple_check_change(expr2, minexpr)
        else:
            fold, change = 'NaN', 0
    return fold, change
        
    
        
        


genes = [];
with open(args.path) as f:
    labels = next(f).strip().split("\t")[1:]
    expr_list = [[] for i in range(len(labels))];
    #print(expr_list)
    for l in f:
        a = l.strip().split("\t")
        genes.append(a)
        for c, el in enumerate(a[1:]):
            expr_list[c].append([float(x) for x in el.split(",")])
        #print(expr_list[0])
        #break
            #res.append( [float(x) for x in el.split(",")] )
        #expr_list.append(res)
        #print(res)
        
#print(len(expr_list[0]))

labels_combinations = list(combinations(labels, 2))
fold_list = [[] for x in labels_combinations]
#print(labels_combinations)
#sys.exit()
for c, (elist1, elist2) in enumerate(combinations(expr_list, 2)):
    for expr1, expr2 in zip(elist1, elist2):
        fold, change = compare(expr1, expr2, args.minfold, args.minexpr)
        score = assign_score(expr1, expr2, fold, change)
        fold_list[c].append((fold, score))
        
 
header = ['gene', 'pair', 'score'] + labels + ["|".join(x) for x in labels_combinations]
print("\t".join(header))
for c, gene in enumerate(genes):
    fold_list_local = [x[c] for x in fold_list]
    score, pos = gene_total_score(fold_list_local)
    label_pair = labels_combinations[pos]
    a = [gene[0], "|".join(label_pair), "%d" % score] + gene[1:] + ["%1.2f" % x[0] if x[0] != 'NaN' else 'NaN' for x in fold_list_local]
    print("\t".join(a))
    
    #print(gene)
    ##print(fold_list_local)
    #print(score, label_pair, pos)
    #print()
    #print()
    
        
        
        #if(change):
            #print(expr1, expr2, fold, change)
    #break;
        



    
        
#size = len(args.path)
#pairs = list(permutations(range(size), 2))
#print(pairs)

#condition_names = [os.path.basename(x).split(".")[0] for x in args.path]
#res_total, stat_total_counts = find_shared_peaks(args.path, args.maxd)
#sys.stderr.write(shared_peaks_stat_to_string(stat_total_counts, size))

#pos2pairs = defaultdict(list);
#for compiled in res_total:
    #for p, peak1, peak2 in compare_multiple_peaks(compiled, pairs, args.minfold, args.mincov):
        #pos2pairs[p].append((peak1, peak2))
    
#for p in pairs:
    #name = "%s_GREATER_%s.gff" % (condition_names[p[0]], condition_names[p[1]])
    #with open(os.path.join(args.outdir, name), 'w') as f:
        #for peak1, peak2 in pos2pairs[p]:
            #f.write(str(pair2interval(peak1, peak2)))
