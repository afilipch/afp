#! /usr/bin/python
'''Annotates the provided genomic regions with the distances to the closest transcription start sites'''

import argparse
import sys
import os
from collections import defaultdict
from bisect import bisect_right, bisect_left

import numpy as np;
from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Annotates the provided genomic regions with the distances to the closest transcription start sites');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genomic regions, gff format");
parser.add_argument('--tss', nargs = '?', required=True, type = str, help = "Path to the TSS annotation file");
parser.add_argument('--stranded', nargs = '?', const=True, default=False, type = bool, help = "If set, the data considered to be stranded");
args = parser.parse_args();


if(os.stat(args.path).st_size != 0):
    ###Connect tss to gene name
    ###################################################################################################
    tss_intervals = BedTool(args.tss);
    gene_starts = defaultdict(list)
    gene_names = defaultdict(list)
    tss = defaultdict(list)
    
    for interval in tss_intervals:
        if(interval[2] == 'TSS'):
            tss[interval.strand].append(interval.start);
        if(interval[2] == 'gene'):
            if(interval.strand == '-'):
                gene_starts[interval.strand].append(interval.end);
            else:
                gene_starts[interval.strand].append(interval.start);
            gene_names[interval.strand].append(interval.attrs['gene_id']);
            
    tss2gene_name = {};
    for strand, tss_starts in tss.items():
        curstarts = list(sorted(gene_starts[strand]))
        curnames = gene_names[strand]
        #print(len(curnames))
        if(strand == '+'):
            for tss_start in tss_starts:
                name = curnames[min(bisect_left(curstarts, tss_start), len(curnames)-1)]
                tss2gene_name[(tss_start, strand)] = name
        else:
            for tss_start in tss_starts:
                name = curnames[max(bisect_left(curstarts, tss_start)-1, 0)]
                tss2gene_name[(tss_start, strand)] = name            
            
    #for el in tss2gene_name.items():
        #print(el)
        
    ###Connect genomic regions to the closest tss;
    ###################################################################################################
    def find_tss(region, tss, strand):
        center = (region.end + region.start)/2
        tss_starts = tss[strand]
        bpos = bisect_left(tss_starts, center)
        tss_up = tss_starts[bpos-1]
        tss_down = tss_starts[bpos]
        if(strand == '+'):
            return (tss_up, center - tss_up), (tss_down, tss_down - center)
        else:
            return (tss_down, tss_down - center), (tss_up, center - tss_up)
    
    
    regions = BedTool(args.path)
    
    
    
    
    
            
else:
    sys.stderr.write("WARNING: the given gff file is empty")







