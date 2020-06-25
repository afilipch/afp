#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Adds motif presense to the all-in table'''


import argparse
import sys
import os




parser = argparse.ArgumentParser(description='Adds motif presense to the all-in table')#, formatter_class = argparse.RawTextHelpFormatter);
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the all in table");
parser.add_argument('--fimo', nargs = '?', required=True, type = str, help = "Path to the fimo output, tsv file")
args = parser.parse_args();


found = [];
with open(args.fimo) as f:
    next(f)
    for l in f:
        a = l.strip().split("\t")
        if(len(a) > 3):
            found.append(a[2])


with open(args.path) as f:
    a = next(f).strip().split("\t")
    a.insert(10, "Motif")
    print("\t".join(a))
    for l in f:
        a = l.strip().split("\t")
        key = "%s:%s:%s" % tuple(a[:3])
        motif = int(key in found)
        #print(motif)
        a.insert(10, str(motif))
        print("\t".join(a))
        #print(">%s:%s:%s\n%s" % (a[0], a[1], a[2], a[9]));
        
    
        
