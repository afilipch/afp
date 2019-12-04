#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts ramified structure of the Prophage database into a single table'''

import argparse
import os
import sys
#import copy
from collections import namedtuple
from glob import glob


from Bio import Entrez
import numpy as np;




parser = argparse.ArgumentParser(description='Converts ramified structure of the Prophage database into a single table');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "UNIX-like pattern path to the prophage features files");
args = parser.parse_args();




### Reading the input

for c, fpath in enumerate(glob(args.path)):
    a = os.path.basename(fpath).strip().split(".")[0].split("_")
    bacteria_id = a[0];
    s, e = int(a[1]), int(a[2])
    if(e>s):
        end = e;
        start = s;
        strand = '+'
    else:
        end = s;
        start = e;
        strand = '-'
        
    num_of_genes = 0;
    integrase = False;
    with open(fpath) as f:
        for l in f:
            fun_ann = l.strip().split("\t")[4]
            if('integrase' == fun_ann or 'Integrase' == fun_ann):
                integrase = True;
                num_of_genes += 1;
                #sys.stderr.write(l);
        num_of_genes += 1;
        
    sys.stdout.write("%s\t%d\t%d\t%d\t%s\t%s\n" % (bacteria_id, start, end, num_of_genes, strand, integrase));
            
    
