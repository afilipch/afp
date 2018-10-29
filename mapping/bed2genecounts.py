#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Assigns expression (TPM, RPKM/FPKM) to the genes based on provided mappings (sorted bed file)'''

import argparse
import os
import sys
from collections import defaultdict, namedtuple
import math

import numpy as np;
import pandas as pd;
import scipy 
import matplotlib.pyplot as plt;
from scipy.stats.stats import pearsonr;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Assigns expression (TPM, RPKM/FPKM) to the genes based on provided mappings (sorted bed file)');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the genomic coverage, bed3 format. First file must be a coverage for plus strand, second for minus");
parser.add_argument('--genes', nargs = '?', required=True, type = str, help = "Path to the genes, gff format");
parser.add_argument('--plot', nargs = '?', type = str, help = "Output destination for the statistics plots");
args = parser.parse_args();

GeneCount = namedtuple('GeneCount', ['genename', 'maxcov', 'variation', 'tpm'])

def assign_coverage(gene, genecov, normfactor):
    if(any(genecov)):
        return GeneCount(genename=gene.name, maxcov=max(genecov), variation = scipy.stats.variation(genecov), tpm = genecov.mean()/normfactor)
    else:
        return GeneCount(genename=gene.name, maxcov=0, variation = 0, tpm = 0)

coverage = {}
coverage['+'] = pd.read_csv(args.path[0], sep="\t" , names = ["chr", "position", "coverage"]).coverage.values
coverage['-'] = pd.read_csv(args.path[1], sep="\t" , names = ["chr", "position", "coverage"]).coverage.values
genes = BedTool(args.genes)

#name2genecount = {};
total_genomic_coverage = sum(sum(x) for x in coverage.values())
#total_gene_coverage = sum([sum(coverage[x.strand][x.start:x.end]) for x in genes])
rpk_sum = 0;
genecounts = [];


for gene in BedTool(args.genes):
    genecov = coverage[gene.strand][gene.start:gene.end];
    rpk_sum += genecov.mean()

rpk_sum = rpk_sum/1000000

print("\t".join(('name', 'tpm', 'maxcov', 'variation')));
for gene in BedTool(args.genes):
    genecov = coverage[gene.strand][gene.start:gene.end];
    gc = assign_coverage(gene, genecov, rpk_sum)
    print("%s\t%.2f\t%d\t%.2f" % (gc.genename, gc.tpm, gc.maxcov, gc.variation))
    genecounts.append(gc)
    
