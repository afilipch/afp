#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Normalizes genomic coverage for bed-pileup file'''

import argparse
import os
import sys
from collections import defaultdict, namedtuple
import math

import numpy as np;
import pandas as pd;


parser = argparse.ArgumentParser(description='Normalizes genomic coverage for bed-pileup file');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the genomic coverage, bed3 format");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Output destination for the statistics plots");
args = parser.parse_args();

total = 0;

for path in args.path:
    total += sum(pd.read_csv(path, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values)

sys.stderr.write("Total coverage: %d\n" % total);    
normfactor = (10**10)/total
sys.stderr.write("Norm-factor: %.2f\n\n" % normfactor);    

for path in args.path:
    outpath = os.path.join(args.outdir, ".".join(os.path.basename(path).split(".")[:-1] + ['normed.bed']))
    sys.stderr.write("Current output: %s\n" % outpath);
    table = pd.read_csv(path, sep="\t", names = ["chr", "position", "coverage"])
    with open(outpath, 'w') as f:
        for row in table.itertuples(index=False):
            f.write("%s\t%d\t%d\n" % (row[0], row[1], row[2]*normfactor))
    
