
'''Adds information from the custom annotation table to the ncbi-like annotated genes'''

import argparse
import sys
import os
from collections import defaultdict
#from bisect import bisect_right, bisect_left

import numpy as np;
from pybedtools import BedTool
from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Adds information from the custom annotation table to the ncbi-like annotated genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genes/transcripts, ncbi-like gff format");
parser.add_argument('--annotation', nargs = '?', required = True, type = str, help = "Path to the custom annotation table, tsv format");
args = parser.parse_args();

ann_dict = {}
cg_dict = {};
with open(args.annotation) as f:
    next(f)
    for l in f:
        a = l.strip().split(";")
        #print(a)
        if(len(a) > 9 and a[1] and a[2]):
            ann_dict[a[1]] = a[2], a[8], a[9]
        if(len(a) > 1 and a[1] and a[0].startswith("cg") and a[1].startswith("NCgl")):
            cg_dict[a[1]] = a[0]
            
#print(len(cg_dict));



for interval in BedTool(args.path):
    name = interval.name.split("-")[1]
    interval.name = name
    interval.attrs['cg'] = cg_dict.get(name, 'None')
    
    lann = ann_dict.get(name, None);
    if(lann):
        genesymbol, annotation, function = lann;
        interval.attrs['genesymbol'] = genesymbol
        interval.attrs['annotation'] = annotation
        interval.attrs['product'] = function
    else:
        interval.attrs['function'] = interval.attrs['product'];
    sys.stdout.write(str(interval))
