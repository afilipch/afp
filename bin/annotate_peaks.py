#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Annotates the discovered peaks'''

import argparse
import sys
import os
from collections import defaultdict
from bisect import bisect_right, bisect_left


import pandas as pd;
import numpy as np;
import matplotlib.pyplot as plt;
from pybedtools import BedTool

from afbio.pybedtools_af import construct_gff_interval
from afbio.numerictools import get_closest

parser = argparse.ArgumentParser(description='Annotates the discovered peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the detected peaks");
parser.add_argument('--coverage', nargs = '?', required=True, type = str, help = "Path to the coverage file");
parser.add_argument('--flen', nargs = '?', default=100, type = int, help = "Length of the peak\'s flanks to be included into analyses");
parser.add_argument('--genes', nargs = '?', required=True, type = str, help = "Path to the gene annotation for the analyzed organizm");

args = parser.parse_args();


###Annotate coverage features of the detected peaks
def annotate_coverage(interval, coverage, flen, normfactor):
    top = int(interval.name)
    start = max((top-flen,0))
    end = top+flen +1
    loccov = coverage[start:end]
    topval = loccov[flen]
    ldrop = min(loccov[:flen])/topval
    rdrop = min(loccov[flen+1:])/topval
    
    return topval/normfactor, ldrop, rdrop, start, end
    
coverage = pd.read_csv(args.coverage, sep="\t" , names = ["chr", "postion", "coverage"]).coverage.values
normfactor = np.mean(coverage);


peaks = BedTool(args.path);
peak2cov = {}
for interval in peaks:
    peak2cov[interval.name] = annotate_coverage(interval, coverage, args.flen, normfactor)
    
    
###Annotate peaks with overlapping genomic features    
genes = BedTool(args.genes)    
rawgenes = peaks.intersect(genes, wo = True)
peak2genes = defaultdict(list);

for el in rawgenes:
    peak2genes[el.name].append(el[9].strip())
    
    
###Get closest genomic starts upstream to the detected peaks
plusstarts = [];
minusstarts = [];
start2genes = defaultdict(list)
for gene in genes:
    if(gene.strand == '+'):
        start = gene.start
        plusstarts.append(start);
    else:
        start = gene.end-1
        minusstarts.append(start);
    start2genes[start].append(gene.name)
minusstarts.sort()
plusstarts.sort()

        
peak2genestarts = {}
for interval in peaks:
    top = int(interval.name)
    cs_plus = plusstarts[min(len(plusstarts)-1, bisect_left(plusstarts, top))];
    d_plus = cs_plus - top
    cs_minus = minusstarts[max(0, bisect_right(minusstarts, top) - 1)];
    d_minus = top - cs_minus
    if(d_plus>d_minus):
        d = d_minus;
        s_strand = '-';
        cs = cs_minus
    else:
        d = d_plus
        s_strand = '+'
        cs = cs_plus
    
    peak2genestarts[interval.name] = start2genes[cs], d, s_strand


for interval in peaks:
    topcoverage, ldrop, rdrop, start, end = peak2cov[interval.name]
    lg = peak2genes[interval.name]
    sg, dg, s_strand = peak2genestarts[interval.name]
    anint = construct_gff_interval(interval.chrom, start, end, 'peak', score=interval.score, strand=interval.strand, source='af_peak_detection', frame='.', attrs=[ ('Name', interval.name), ('topcoverage', "%1.3f" % topcoverage), ('ldrop', "%1.3f" % ldrop), ('rdrop', "%1.3f" % rdrop), ('genes', ",".join(lg)), ('start_genes', ",".join(sg)), ('start_gene_distance', dg), ('start_gene_strand', s_strand) ])
    sys.stdout.write(str(anint))
#construct_gff_interval(chrom, start, stop, feature, score='0', strand='.', source='un', frame='.', attrs=[]);




