#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Checks binding enrichment inside the given locus for multiple binding factors, creates tables'''

import argparse
import os
import sys
from os import listdir
from os.path import isfile
import numpy as np;
from collections import defaultdict

from pybedtools import BedTool, Interval





parser = argparse.ArgumentParser(description='Checks binding enrichment inside the given locus for multiple binding factors, creates tables');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with binding peaks. \'peaks\' folder in chipchap project folder");
parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Path to experiments annotation table");
parser.add_argument('--start', nargs = '?', required=True, type = int, help = "Start of the region");
parser.add_argument('--end', nargs = '?', required=True, type = int, help = "End of the region");
parser.add_argument('--length', nargs = '?', required=True, type = int, help = "Length of the genome");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();

region = BedTool([Interval('chr1', args.start, args.end, strand = '+', score = '0', name = 'region')])

exp2protein = {};
with open(args.table) as f:
    next(f);
    for l in f:
        a = l.strip().split("\t")
        exp2protein[".".join(a[1].split(".")[:-1])] = a[3]
#print(exp2protein)
    
peakfiles = [os.path.join(args.path, f) for f in listdir(args.path) if isfile(os.path.join(args.path, f)) and 'annotated' in f]
protein2exp = defaultdict(list)
for pf in peakfiles:
    name = ".".join( os.path.basename(pf).split(".")[:-2] )
    protein2exp[exp2protein[name]].append(pf)    



def get_enrichment(peakfile, region, length):
    peaks = BedTool(peakfile)
    ilen = len(region[0])
    olen = length - ilen
    outside = [float(x.attrs['topcoverage']) for x in peaks.intersect(region, f = 0.5, v = True)]
    inside = [float(x.attrs['topcoverage']) for x in peaks.intersect(region, f = 0.5, u = True)]
    
    if(len(outside)):
        normed_count = (len(inside)/ilen)/(len(outside)/olen)
        normed_sum = (sum(inside)/ilen)/(sum(outside)/olen)
        overrepresented = int(normed_count > 2 or normed_sum > 2)
    else:
        normed_count = float('nan')
        normed_sum = float('nan')
        overrepresented = -1
    return "%s\t%d\t%d\t%1.1f\t%1.1f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%d\n"  % (os.path.basename(peakfile).split(".")[0] , len(inside), len(outside), sum(inside), sum(outside), np.mean(inside), np.mean(outside), np.median(inside), np.median(outside), normed_count, normed_sum, overrepresented)


for protein, peaklist in protein2exp.items():
    with open(os.path.join(args.outdir, "%s_prophage_enrichment.tsv" % protein), 'w') as f:
        f.write("experiment\tnumber peaks inside\tnumber peaks outside\ttotal peak coverage inside\ttotal peak coverage outside\tmean peak coverage inside\tmean peak coverage outside\tmedian peak coverage inside\tmedian peak coverage outside\tpeak density ration\tpeak total coverage ratio\toverrepresented\n")
        for peakfile in peaklist:
            f.write(get_enrichment(peakfile, region, args.length));





