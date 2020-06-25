#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Extracts sequences from table for FIMO analyses'''


import argparse
import sys
import os




parser = argparse.ArgumentParser(description='Extracts sequences from table for FIMO analyses')#, formatter_class = argparse.RawTextHelpFormatter);
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the all in table");
args = parser.parse_args();


with open(args.path) as f:
    next(f)
    for l in f:
        a = l.strip().split("\t")
        print(">%s:%s:%s\n%s" % (a[0], a[1], a[2], a[9]));
        
    
        
