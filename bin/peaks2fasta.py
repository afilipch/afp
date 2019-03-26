#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts annotated peaks in gff to fasta format compatible with meme'''

import argparse
import os
import sys
import math
#import scipy
import numpy as np;
import pandas as pd;
from pybedtools import BedTool
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Converts peaks in bed/gff to fasta format compatible with meme');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the peaks");
parser.add_argument('--genome', nargs = '?', required = True, type = str, help = "Path to the genome, fasta format");
parser.add_argument('--weighted', nargs = '?', default=False, const=True, type = bool, help = "If set, coverage vs convolution is plotted");
args = parser.parse_args();

genome = next(SeqIO.parse(args.genome, 'fasta')).seq
if(args.weighted):
    scores = [math.log(float(x.attrs['topcoverage'])+1) for x in BedTool(args.path)]
    norm = max(scores);
    scores = [x/norm for x in scores];
    header = '>WEIGHTS %s' % ' '.join(['%1.3f' % x for x in scores])
    print(header);
    
              
for interval in BedTool(args.path):
    seq = genome[interval.start: interval.end].upper()
    sid = interval.name
    print(">%s\n%s" % (sid, seq))
