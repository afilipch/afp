#! /usr/local/anaconda3/bin/python
'''Generates artificial peaks'''

import argparse
import os
import sys
import numpy as np;
import random

from Bio import SeqIO
from pybedtools import BedTool

from afbio.generators import generator_doublesam
from afbio.numerictools import select_by_probability
import matplotlib.pyplot as plt;

parser = argparse.ArgumentParser(description='Generates artificial peaks');

parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to genome, fasta file");

parser.add_argument('--numpeaks', nargs = '?', default=400, type = int, help = "Number of reads to generate");
parser.add_argument('--doublets', nargs = '?', default=False, const = True, type = bool, help = "Number of samples to generate");
parser.add_argument('--doublet_ratio', nargs = '?', default=0.5, type = float, help = "Ratio between bigger and smaller peaks");
args = parser.parse_args();


DRANGE = list(range(10, 100));

def score_func(x):
    return 2**(-x)

score_range = np.linspace(0, 10, 1001)
score_probs = [score_func(x) for x in score_range]
scores = [3 + x for x in random.choices(score_range, score_probs, k=args.numpeaks)]

genome = dict([(x.name, len(x)) for x in SeqIO.parse(args.genome, 'fasta')]);
chrnames, chrsizes = list(zip(*genome.items()))
 

for score in scores:
    chrom = select_by_probability(chrnames, chrsizes);
    size = genome[chrom];
    pos = random.randint(1000, size - 1000);
    print("%s\t%d\t%d\t%d\t%1.3f\t+" % (chrom, pos, pos+1, pos, score))
    if(args.doublets):
        dpos = pos + random.choice((1,-1))*random.choice(DRANGE)
        dscore = score*args.doublet_ratio
        print("%s\t%d\t%d\t%d\t%1.3f\t+" % (chrom, dpos, dpos+1, dpos, dscore))
