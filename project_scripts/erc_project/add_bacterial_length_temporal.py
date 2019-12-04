#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts ramified structure of the Prophage database into a single table'''

import argparse
import os
import sys
#import copy
from collections import namedtuple
from glob import glob


from Bio import SeqIO




parser = argparse.ArgumentParser(description='adds an information on the length of bacteria');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "UNIX-like pattern path to the bacteria contigs (seeds) files");
parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Path to the bacteria metainformation table (prophage_metainfo/bacteria2taxonomy.tsv)");
args = parser.parse_args();




### Reading the input

chrom2length = {};
nums = set()
for c, fpath in enumerate(glob(args.path)):
    fasta_path = os.path.join(fpath, 'contigs')
    num = 0;
    for seqrecord in SeqIO.parse(fasta_path, 'fasta'):
        chrom2length[seqrecord.name] = len(seqrecord);
        num+=1;
    nums.add(num)
    if(c and c % 2000 == 0):
        sys.stderr.write("%d files processed\n" % c);
        


with open(args.table) as f:
    for l in f:
        a = l.strip().split("\t");
        bactera_id = a[0]
        length = chrom2length.get(bactera_id, None)
        if(length):
            a.append(str(length));
            print("\t".join(a))
        else:
            sys.stderr.write("For bacteria %s infromation is missing\n" % bactera_id)
            

sys.stderr.write("\n%s\n\n" % nums)