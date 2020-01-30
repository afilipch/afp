#! /usr/bin/python
'''Converts annotated transcripts (with alternative 5'UTRS and 3'UTRS) into bed12 format compatible with UCSC'''

import argparse
import sys
import os
from collections import defaultdict, Counter

import numpy as np;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.numerictools import CDF


parser = argparse.ArgumentParser(description='Converts transcripts (with alternative 5\'UTRS) into bed12 format compatible with UCSC');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated transcripts file, gff format");
parser.add_argument('--trackopts', nargs = '?', default='', type = str, help = "Header for the (ucsc) track")
args = parser.parse_args();


if(args.trackopts):
    print(args.trackopts);




def block2bed12(block):
    for c, gene in enumerate(block, start=1):
        cds_start, cds_stop = gene.attrs['cds'].split(":")
        print("\t".join(('chr1', str(gene.start), str(gene.stop), "%s_%d" % (gene.name, c), gene.score, gene.strand, cds_start, cds_stop, '0', '1', str(len(gene)), '0')))
    
    

genes =  defaultdict(list)
for x in BedTool(args.path):
    genes[x.name].append(x)




for block in genes.values():
    block2bed12(block)
