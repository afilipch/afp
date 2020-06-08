#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Checks the recovery and false discovery rate of Peakaboo'''

import argparse
import os
import sys
from itertools import combinations;



import numpy as np;
from scipy.stats import pearsonr, normaltest, variation
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;



parser = argparse.ArgumentParser(description='Checks the recovery and false discovery rate of Peakaboo');
parser.add_argument('--original', nargs = '?', required=True, type = str, help = "Path to the original peaks, bed file");
parser.add_argument('--detected', nargs = '?', required=True, type = str, help = "Path to the detected annotated peaks, gff format");
#parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
#parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();


def find_closest(intervals_1, intervals_2):
    res = []
    start = 0;
    for i1 in intervals_1:
        curint = intervals_2[start]
        curd = abs(i1[0] - curint[0])
        for p, i2 in enumerate(intervals_2[start+1:]):
            d = abs(i1[0] - i2[0])
            if(d<curd):
                curd = d;
                curint = i2
            else:
                res.append((i1, curint, curd))
                start += p
                break;
        else:
            res.append((i1, curint, curd))
            start += p
    return(res)
                
        


original = [(int(x.name), float(x.score)) for x in BedTool(args.original)]
detected = [(int(x.name), float(x.attrs['topcoverage'])) for x in BedTool(args.detected)]
original.sort(key = lambda x: x[0])
detected.sort(key = lambda x: x[0])



recovery = find_closest(original, detected);
print(pearsonr([x[0][1] for x in recovery], [x[1][1] for x in recovery]))


print(recovery[0])

#for el in recovery[:100]:
    #print(el)

