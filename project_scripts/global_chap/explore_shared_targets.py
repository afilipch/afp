#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores gene networks for mutliple binding proteins'''

import argparse
import os
import sys
from os import listdir
from os.path import isfile
from collections import defaultdict

import numpy as np;
from scipy.stats import percentileofscore
from pybedtools import BedTool, Interval
from itertools import product

from afbio.pybedtools_af import read_comments




parser = argparse.ArgumentParser(description='Explores shared targets among different dna-binding proteins');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the merged binding regions, gff files");
#parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to NCBI annotation table");
#parser.add_argument('--maxsd', nargs = '?', default = 350, type = int, help = "Maximal allowed distance to the start gene");
#parser.add_argument('--strict', nargs = '?', default = False, const=True, type = bool, help = "If set, strict requirments for a gene to be counted as controlled by an RBP are applied");
args = parser.parse_args();



position2name = dict([ (x[0], x[1].split("_")[0]) for x in enumerate(read_comments(args.path)[0][2:].split(",")) ])
order = list(set(position2name.values()));
order.sort()


intensities = [ [float(y) for y in x.attrs['topcoverage'].split(",")] for x in BedTool(args.path)]

shared_list = []
for interval in BedTool(args.path):
    lints = [float(x) for x in interval.attrs['topcoverage'].split(",")]
    names = [position2name[x[0]] for x in enumerate(lints) if x[1]]
    shared_list.append(names)
    if(len(set(names)) > 1):
        sys.stderr.write(str(interval))

for n1 in order:
    intersections = []
    for n2 in order:
        print("%s\t%s\t%d\n" % (n1, n2, len([x for x in shared_list if n1 in x and n2 in x])) )
#for el in shared_list:
    #print(el)
 
#print(len(shared_list))
#print( [x for x in shared_list if len(list(set(x)))>1] )
    
##for el in shared_list:
    ##print(el);




    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
