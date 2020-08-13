#! /usr/bin/python
'''Converts transcripts (with alternative 5'UTRS) into bed12 format compatible with UCSC'''

import argparse
import sys
import os
from collections import defaultdict, Counter

import numpy as np;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.numerictools import CDF


parser = argparse.ArgumentParser(description='Converts transcripts (with alternative 5\'UTRS) into bed12 format compatible with UCSC');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the transcripts file, gff format");
#parser.add_argument('--cds', nargs = '?', default='', type = str, help = "Path to the CDS, gff bed format")
parser.add_argument('--trackopts', nargs = '?', default='', type = str, help = "Header for the (ucsc) track")
args = parser.parse_args();

if(args.trackopts):
    print(args.trackopts);
    
    
chr2region = defaultdict(list)
for region in BedTool(args.path):
    chr2region[region.chrom].append(region)
    
chr2region = [x[1] for x in sorted(chr2region.items(), key = lambda y: y[0])]

for chr_count, region_list in enumerate(chr2region, start = 1):
    for region in region_list:
        print("\t".join(('chr%d' % chr_count, str(region.start), str(region.end), region.name, region.score, region.strand, str(region.start), str(region.stop), '0', '1', str(len(region)), '0')))
    


#if(args.trackopts):
    #print(args.trackopts);


#def block2bed12(name, block, cds_dict):
    #cds = cds_dict.get(name, None)
    #for n, region in enumerate(block, start=1):
        #if(not cds):
            #cds = region
        #print("\t".join(('chr1', str(region.start), str(region.end), "%s_%d" % (region.name, n), region.score, region.strand, str(cds.start), str(cds.stop), '0', '1', str(len(region)), '0')))
    
    

#genes =  defaultdict(list)
#for x in BedTool(args.path):
    #genes[x.name].append(x)

#if(args.cds):
    #cds_dict = dict([ (x.name, x) for x in  BedTool(args.cds) ])
#else:
    #cds_dict = {}

#for name, block in genes.items():
    #block2bed12(name, block, cds_dict)
