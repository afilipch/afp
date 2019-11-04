#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Adds divergent targets to the hrra all in table'''

import argparse
import os
import sys

import numpy as np
from pybedtools import BedTool, Interval
from collections import defaultdict, Counter


parser = argparse.ArgumentParser(description='Adds divergent targets to the hrra all in table');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the hrra all-in table");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the custom annotation with gene names");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts with tss");
parser.add_argument('--cds', nargs = '?', required=True, type = str, help = "Path to the CDS");
parser.add_argument('--diff', nargs = 3, required=True, type = str, help = "Path to the file with genes\' differential expression");
args = parser.parse_args()


def get_divergent(peak, starts, stops, strand, gene2cds, maxd=300, inside = 100):
    if(strand=='-'):
        distances = [(x[0], x[1] - peak) for x in starts]
    else:
        distances = [(x[0], peak - x[1]) for x in stops]
    distances = [x for x in distances if x[1]>-1*inside]
    gene_name, mindistance = min(distances, key = lambda x: abs(x[1]))
    if(mindistance<=maxd):
        cds_start, cds_stop = gene2cds[gene_name]
        if(strand=='-'):
            atg = cds_start - peak;
        else:
            atg = peak - cds_stop
        return gene_name, mindistance, atg
    else:
        return None
        


tr_list = BedTool(args.transcripts)
starts_plus = [ (x.attrs['Parent'].split("-")[1], x.start) for x in tr_list if x.strand == '+']
stops_minus = [ (x.attrs['Parent'].split("-")[1], x.stop) for x in tr_list if x.strand == '-']
gene2strand = dict( [ (x.attrs['Parent'].split("-")[1], x.strand) for x in tr_list] )

cds_list = BedTool(args.cds)
gene2cds = dict( [ (x.attrs['Parent'].split("-")[1], (x.start, x.stop) ) for x in cds_list] )

gene2annotation = {};
gene2genesymbol = {};
with open(args.annotation) as f:
    next(f);
    for l in f:
        a = l.strip().split("\t")
        
        if(len(a) > 8 and a[1]):
            if(a[2]):
                gene2genesymbol[a[1]] = a[2]
            else:
                gene2genesymbol[a[1]] = a[1]
            gene2annotation[a[1]] = a[7], a[8]
                


gene2diff = defaultdict(list)
for path in args.diff:
    with open(path) as f:
        next(f)
        for l in f:
            a = l.strip().split("\t");
            name = a[0].split("-")[1]
            wt = np.mean([float(x) for x in a[1].split(";")])
            ko = np.mean([float(x) for x in a[2].split(";")])
            log_ko_wt = float(a[3]);
            
            gene2diff[name].append(("%1.1f" % wt,"%1.1f" %  ko,"%1.3f" %  log_ko_wt))



with open(args.path) as f:
    a = next(f).strip().split("\t")
    #print(a[14:])
    a.append("divergent");
    print("\t".join(a))
    for l in f:
        a = l.strip().split("\t")
        a.append("no");
        print("\t".join(a))
        start, end = int(a[1]), int(a[2])
        peak = (start+end)//2
        strand = gene2strand[a[3]]
        res = get_divergent(peak, starts_plus, stops_minus, strand, gene2cds ,maxd=300, inside = 100)
        if(res):
            gene_name, mindistance, atg = res;
            function, annotation = gene2annotation[gene_name]
            genesymbol = gene2genesymbol[gene_name]
            a[3] = gene_name
            a[4] = genesymbol
            a[5] = str(atg)
            a[6] = str(mindistance)
            a[-3] = function
            a[-2] = annotation
            a[-1] = "yes"
            
            diff = gene2diff[gene_name]
            a[14] = diff[0][0]
            a[15] = diff[1][0]
            a[16] = diff[2][0]
            
            a[17] = diff[0][1]
            a[18] = diff[1][1]
            a[19] = diff[2][1]
            
            a[20] = diff[0][2]
            a[21] = diff[1][2]
            a[22] = diff[2][2]
            print("\t".join(a))
            
            

