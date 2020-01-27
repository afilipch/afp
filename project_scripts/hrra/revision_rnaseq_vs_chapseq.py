#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Correlates differential expression of the genes with the peak intensity on their promoters [revision version]'''



import argparse
import os
import sys
#import copy
from collections import defaultdict, Counter
from itertools import product;

import numpy as np;
from scipy.stats import sem
#import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Correlates differential expression of the genes with the peak intensity on their promoters [revision version]');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the all-in table");
parser.add_argument('--diff', nargs = '?', required=True, type = str, help = "Path to the file with genes\' differential expression");
parser.add_argument('--distance', nargs = '?', default=400, type = int, help = "Maximum allowed distance for the peak to the closest peak");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Path to the plot");
args = parser.parse_args();



def discriminate(gene2diff, bound_set, labels, start):
    bound, unbound = defaultdict(list), defaultdict(list);
    for geneid, currdiff in gene2diff.items():
        for ldiff, label in zip(currdiff[start:], labels[start:]):
            if(sum(ldiff[:2]) > 10):
                change = abs(ldiff[2]);
                if(change < 100):
                    if(geneid in bound_set):
                        bound[label].append(abs(ldiff[2]))
                    else:
                        unbound[label].append(abs(ldiff[2]))
    return bound, unbound;
                    

bound_genes = set();
with open(args.path) as f:
    next(f)
    for l in f:
        a = l.strip().split("\t")
        geneid = a[3]
        tss = float(a[6]);
        bound_genes.add(geneid);
            
            

gene2diff = {}
with open(args.diff) as f:
    next(f)
    for l in f:
        a = l.strip().split("\t");
        name = a[0].split("-")[1]
        ko = [float(x) for x in a[1].split(";")]
        wt = [float(x) for x in a[2].split(";")]
        if_expr = sum(ko) > 20 or sum(wt) > 20
        if_fold = abs(float(a[4]))>1/3
        
        log_ko_wt = float(a[3]);
        if_diff = bool(int(a[5]))
        gene2diff[name] = log_ko_wt, if_diff and if_fold#and if_expr and if_fold
       
       
       
       
#for k, v in gene2diff.items():
    #print(v)
    
#print(len(set(bound_genes)))

total_true = len([x for x in gene2diff.values() if x[1]])
bound_true = len([x for x in gene2diff.items() if x[1][1] and x[0] in bound_genes])
print("total true: %d\nbound true: %d\n" % (total_true, bound_true))


bound_expression = [gene2diff[x] for x in bound_genes];
total_expression = [x[0] for x in bound_expression]
true_expression = [x[0] for x in bound_expression if x[1]]

total_increase = len([x for x in total_expression if x>0])
total_decrease = len([x for x in total_expression if x<0])
true_increase = len([x for x in true_expression if x>0])
true_decrease = len([x for x in true_expression if x<0])

print("true increase: %d\ntrue decrease: %d\ntotal increase: %d\ntotal decrease: %d\n" % (true_increase, true_decrease, total_increase, total_decrease))


#for geneid in bound_genes:
    ##print(geneid, gene2diff[geneid])
    #print(gene2diff[geneid][1])
    
    

def get_stat_paired(pair):

    if(pair[0].is_reverse):
        strand = '+'
    else:
        strand = '-'
    
    chrom = pair[0].reference_name
    start = min([x.reference_start for x in pair]);
    end = max([x.reference_end for x in pair]);
    score = sum([float(x.get_tag('AS')) for x in pair])
            
    return chrom, start, end, score, strand 


def get_stat_single(s):
    if(s.is_reverse):
        strand = '+'
    else:
        strand = '-'
    
    chrom = s.reference_name
    start = s.reference_start 
    end = s.reference_end
    score = float(s.get_tag('AS'))
            
    return chrom, start, end, score, strand 
    



def readmappings2stat(readmappings, func, ambiguous):
    if(not readmappings):
        return []
    if(not ambiguous and len(readmappings)>1):
            return []
    return [func(x) for x in readmappings]




    





samfile = pysam.AlignmentFile(args.path)
for readmappings in generator_mappings(samfile):
    stat_list = readmappings2stat(readmappings, get_stat, args.ambiguous)


