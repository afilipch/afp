#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores motifs without peaks'''

import argparse
import os
import sys
from collections import defaultdict
import copy
from pathlib import Path


import numpy as np;
from pybedtools import BedTool
from itertools import combinations, permutations
import matplotlib.pyplot as plt;
from afbio.pybedtools_af import read_comments
from afbio.peaks import find_closest_peak_unstranded




parser = argparse.ArgumentParser(description='Explores motifs without peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated (merged) peaks, gff format");
parser.add_argument('--fimo', nargs = '?', required = True, type = str, help = "Path to the fimo output, gff format");
parser.add_argument('--mind', nargs = '?', default = 50, type = int, help = "Minimal allowed distance from motif to the closest peaks, to be counted as orphans");
args = parser.parse_args()



#comments = read_comments(args.path)[0]
#labels = comments.split(' ')[1].split(",");
fimo = BedTool(args.fimo)
peaks = BedTool(args.path)
orphans = []

for feature, peak, distance in find_closest_peak_unstranded(peaks, fimo):
    if(abs(distance)>args.mind):
        orphans.append(feature);
        
print(orphans[2])


    


