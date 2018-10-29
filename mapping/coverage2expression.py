#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Evaluates gene expression values based on coverage'''

import argparse
import os
import sys
from collections import defaultdict, namedtuple
import math
from itertools import combinations

import numpy as np;
import pandas as pd;
import scipy 
import matplotlib.pyplot as plt;
from scipy.stats.stats import pearsonr;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Evaluates gene expression values based on coverage');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "Path to the sam file. For example: coverage2expression s1.plus.bed s1.minus.bed");
parser.add_argument('--genes', nargs = '?', required=True, type = str, help = "Path to the annotated genes, gff format");

args = parser.parse_args();

genecounts = [];

#GeneCount = namedtuple('Gene', ['name', 'counts1', 'counts2', 'fold', 'variation1', 'variation2'])
strand2pos = {'+': 0, '-': 1};

def get_coverage(path):
    return pd.read_csv(path, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values

coverages = [get_coverage(x) for x in args.path]



def gene2counts(gene, coverages):
    if(gene.strand == '+'):
        scov = coverages[0]
    else:
        scov = coverages[1]
        
    return np.mean(scov[gene.start:gene.end])

 
genes = list(BedTool(args.genes))
rawcounts = [gene2counts(x, coverages) for x in genes]

norma = 1000000/sum(rawcounts);
tpms = [x*norma for x in rawcounts];

for gene, tpm in zip(genes, tpms):
    print("%s\t%1.2f" % (gene.name, tpm))

    

        

    
    
    





