#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores gene networks for mutliple binding proteins'''

import argparse
import os
import sys
from os import listdir
from os.path import isfile
from collections import defaultdict

import numpy as np;
from scipy.stats import percentileofscore
from pybedtools import BedTool, Interval

from afbio.sequencetools import get_sub_lists




parser = argparse.ArgumentParser(description='Explores gene networks for mutliple binding proteins');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with binding peaks. \'peaks\' folder in chipchap project folder");
parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Path to experiments annotation table");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to NCBI annotation table");
parser.add_argument('--maxsd', nargs = '?', default = 350, type = int, help = "Maximal allowed distance to the start gene");
parser.add_argument('--strict', nargs = '?', default = False, const=True, type = bool, help = "If set, strict requirments for a gene to be counted as controlled by an RBP are applied");
args = parser.parse_args();

id2gene = dict( [ (x.attrs['ID'], x.attrs['Name']) for x in BedTool(args.annotation) if x[2]=='gene'] )



exp2protein = {};
with open(args.table) as f:
    next(f);
    for l in f:
        a = l.strip().split("\t")
        exp2protein[".".join(a[1].split(".")[:-1])] = a[3]
#print(exp2protein)
    
peakfiles = [os.path.join(args.path, f) for f in listdir(args.path) if isfile(os.path.join(args.path, f)) and 'annotated' in f]
for pf in peakfiles:
    name = ".".join( os.path.basename(pf).split(".")[:-2] )
    exp2protein[pf] = exp2protein[name]


def peakfile2genes(peakfile, maxsd):
    peaks = BedTool(peakfile);
    scores = [float(x.attrs['topcoverage']) for x in peaks]
    gene2scores = defaultdict(float)
    for peak in peaks:
        if(int(peak.attrs['start_gene_distance']) <= maxsd):
            gene2scores[id2gene[peak.attrs['start_gene']]] += float(peak.attrs['topcoverage'])
            
    return dict([ (x[0], percentileofscore(scores, x[1])) for x in gene2scores.items()])


rbp2genesets = defaultdict(list)
for peakfile in peakfiles:
    rbp2genesets[exp2protein[peakfile]].append(peakfile2genes(peakfile, args.maxsd))
    
def collapse_genes(genesets):
    gene2scores = defaultdict(list);
    allgenes = set()
    for d in genesets:
        allgenes.update(d.keys());
    for gene in allgenes:
        for d in genesets:
            gene2scores[gene].append(d.get(gene, 0));
            
    return gene2scores


def strict_filter_genes(gene2scores):
    genes = []
    for gene, scores in gene2scores.items():
        nonzero = [x for x in scores if x];
        if(np.mean(nonzero) >= 50 and (len(nonzero)/len(scores) >= 0.5 or len(nonzero)>=3)):
            genes.append(gene)
    return set(genes)

def soft_filter_genes(gene2scores):
    genes = []
    for gene, scores in gene2scores.items():
        nonzero = [x for x in scores if x];
        if(np.mean(nonzero) >= 25 and (len(nonzero)/len(scores) >= 0.25 or len(nonzero)>=2)):
            genes.append(gene)
    return set(genes)

if(args.strict):
    filter_genes = strict_filter_genes
else:
    filter_genes = soft_filter_genes
    
rbp2genes = {};
for rbp, genesets in rbp2genesets.items():
    if(len(genesets) >= 4):
        
        gene2scores = collapse_genes(genesets)
        rbp2genes[rbp] = filter_genes(gene2scores);
    
for rbpset in get_sub_lists(rbp2genes.keys(), 2):
    set_genes = [rbp2genes[x] for x in rbpset];
    intersection_genes = set.intersection(*set_genes);
    #print(intersection_genes)
    print("%s\t%s" % ( ",".join(list(rbpset)), ",".join(list(intersection_genes)) ) );
#for kv in collapse_genes(genesets).items():
    #print(kv)
    


    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
