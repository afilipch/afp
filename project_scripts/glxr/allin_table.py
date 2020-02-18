#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates table with all the information regarding GLXR camp project'''

import argparse
import os
import sys
from collections import defaultdict, Counter


from Bio import SeqIO

import numpy as np;
from pybedtools import BedTool, Interval

from afbio.pybedtools_af import read_comments




parser = argparse.ArgumentParser(description='Creates table with all the information regarding GLXR camp project');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated (with antisense addtion) consensus regions");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file")
parser.add_argument('--mincov', nargs = '?', default=1, type = float, help = "Minimum allowed coverage")
args = parser.parse_args();

FLANK = 30;
SAMPLE_REPLICATES = [0,3,5,7,9]

sample_names = read_comments(args.path)[0][1:].strip().split(",")
header = ["Chrom", "Start", "Stop", "Gene ID", "Gene symbol", "Distance ATG", "Distance to TSS", "Predicted Function", "Annotation", "Sequence"] + ["TC %s" % x for x in sample_names] + ["AC %s" % x for x in sample_names] + ["AS Gene ID", "AS Gene symbol", "AS Distance ATG", "AS Distance to TSS", "AS Predicted Function", "AS Annotation"]


def interval2seq(interval, reference):
    start = int(interval.name) - FLANK
    stop = int(interval.name) + FLANK +1
    if(interval.strand == '+'):
        return str(reference[interval.chrom][start:stop].seq.upper())
    elif(interval.strand == '-'):
        return str(reference[interval.chrom][start:stop].seq.reverse_complement().upper())


def check_interval(interval, mincov):
    repl_list = [float(x) for x in interval.attrs['topcoverage'].split(",")]
    repl_list = [repl_list[x[0]:x[1]] for x in zip(SAMPLE_REPLICATES, SAMPLE_REPLICATES[1:])]
    if(any([all([y>mincov for y in x]) for x in repl_list])):
        return max([sum(x) for x in repl_list]);
    else:
        return 0;
    
    
def output_interval(interval, antisense, genome):
    a = [interval.chrom, str(interval.start), str(interval.stop), interval.attrs['cg'],  interval.attrs['genesymbol'], interval.attrs['atg'], interval.attrs['tss'], interval.attrs['function'], interval.attrs['annotation'], interval2seq(interval, genome)] + interval.attrs['topcoverage'].split(",") + interval.attrs['area_coverage'].split(",") 
    if(antisense):
        a.extend([antisense.attrs['cg'],  antisense.attrs['genesymbol'], antisense.attrs['atg'], antisense.attrs['tss'], antisense.attrs['function'], antisense.attrs['annotation']])
    else:
        a.extend(["None"]*6)
        
    print("\t".join(a))
    
    
genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))


antisense_dict = {};
passed_intervals = [];

for interval in BedTool(args.path):
    if(interval.attrs['anti'] == '0'):
        score = check_interval(interval, args.mincov)
        if(score):
            passed_intervals.append((interval, score))
    else:
        antisense_dict[(interval.chrom, interval.name)] = interval;


print("\t".join(header))
for interval, _ in sorted(passed_intervals, key = lambda x: x[1], reverse = True):
    antisense = antisense_dict.get( (interval.chrom, interval.name) )
    output_interval(interval, antisense, genome);
    

            
            
