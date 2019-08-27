#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Assigns expression TPM to the genes based on provided mappings (sorted bed file)'''

import argparse
import os
import sys
from collections import defaultdict
import math
import numpy as np;
from scipy.stats import variation 
import matplotlib.pyplot as plt;
#from scipy.stats.stats import pearsonr;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Assigns expression TPM to the genes based on provided mappings (sorted bed file)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the mapped reads, sorted bed format");
parser.add_argument('--genes', nargs = '?', required=True, type = str, help = "Path to the genes, sorted gff format");
parser.add_argument('--stranded', nargs = '?', default=False, const = True, type = str, help = "If set the RNA-seq data are supposed to be stranded");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Output destination for the statistics plots");
args = parser.parse_args();

def normalize(geneid2tval):
    norma = sum(geneid2tval.values())
    return dict([ (x[0], x[1]*1000000/norma) for x in geneid2tval.items() ])

def get_geneid(intersection):
    attrs = dict( [x.split("=") for x in intersection[-2].split(";")])
    return attrs["ID"]

mapped_reads = BedTool(args.path)
genes = BedTool(args.genes);
geneid2tval = {};

if(args.stranded):
    mapped_to_genes = 0
    for cov in genes.coverage(b=mapped_reads, F=0.51, s=True, sorted=True):
        geneid2tval[cov.name] = int(cov[9])/len(cov);
        mapped_to_genes += int(cov[9])
        
    sys.stderr.write("\nTotal reads: %d\nReads mapped to genes: %d\nFraction mapped %1.2f\n\n" % (len(mapped_reads), mapped_to_genes, mapped_to_genes/len(mapped_reads)))
        
        
    
else:
    geneid2count = defaultdict(float);
    curname = ''
    cur_geneids = [];
    shared_reads = 0;
    for interval in mapped_reads.intersect(b = genes, wo = True, f=0.51, sorted=True):
        #print(curname);
        #print(interval.name)
        if(interval.name == curname):
            cur_geneids.append(get_geneid(interval))
        else:
            curname = interval.name;
            if(curname):
                for geneid in cur_geneids:
                    geneid2count[geneid] += 1/len(cur_geneids)
                if(len(cur_geneids) > 1):
                    shared_reads += 1;
            cur_geneids = [get_geneid(interval)];
            
    for gene in genes:
        geneid2tval[gene.name] = geneid2count[gene.name]/len(gene);
        
    sys.stderr.write("\nTotal reads: %d\nReads mapped ambiguosly: %d\nFraction of ambiguosly mapped %1.2f\n\n" % (len(mapped_reads), shared_reads, shared_reads/len(mapped_reads)))
        
        
        
geneid2tpm = normalize(geneid2tval)
for gene in genes:
    gene.attrs['tpm'] = "%1.2f" % geneid2tpm[gene.name]
    sys.stdout.write(str(gene))
    
    
sys.stderr.write("\nSample: %s\nTotal number of genes: %d\nNumber of expressed genes: %d\n\n" % (os.path.basename(args.path).split(".")[0], len(genes), len([x for x in geneid2tval.values() if x])) )
        
        
        
        
        
        
        
        
        
        
    
