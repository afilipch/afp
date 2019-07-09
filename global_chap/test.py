#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores gene networks for mutliple binding proteins'''

import argparse
import os
import sys
from os import listdir
from os.path import isfile
import numpy as np;

from pybedtools import BedTool, Interval





parser = argparse.ArgumentParser(description='Explores gene networks for mutliple binding proteins');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with binding peaks. \'peaks\' folder in chipchap project folder");
parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Experiments annotation table");
#parser.add_argument('--end', nargs = '?', required=True, type = int, help = "End of the region");
#parser.add_argument('--length', nargs = '?', required=True, type = int, help = "Length of the genome");
args = parser.parse_args();

exp2protein = {};
with open(args.table) as f:
    next(f);
    for l in f:
        a = l.strip().split("\t")
        exp2protein[".".join(a[1].split(".")[:-1])] = a[3]
#print(exp2protein)
    
peakfiles = [os.path.join(args.path, f) for f in listdir(args.path) if isfile(os.path.join(args.path, f)) and 'annotated' in f]

for pf in peakfiles:
    name = ".".join( os.path.basename(pf).split(".")[:-2] )
    if name not in exp2protein:
        print(name)
    #print(pf, exp2protein.get(name, None));
    
#for k in  exp2protein:
    #print(k);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
