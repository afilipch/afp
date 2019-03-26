#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Provides basic description of the input fasta file'''

import argparse
import os
import sys


from Bio import SeqIO
from collections import Counter

parser = argparse.ArgumentParser(description='Provides basic description of the input fasta file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the fasta file");
#parser.add_argument('--genes', nargs = '?', required=True, type = str, help = "Path to the genes, sorted gff format");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Output destination for the statistics plots");
args = parser.parse_args();

NUCLS = 'ACTGN'
def get_basic_stat(seqrecord):
    length = len(seqrecord);
    nucl_counter = Counter(str(seqrecord.seq).upper())
    nucl_freq = ["%.3f" % (nucl_counter.get(x, 0)/float(length)) for x in NUCLS];
    return length, nucl_freq, nucl_counter
        

header = ["Name", "Length", "A", "C", "T", "G", "N"]
print("\t".join(header))
total_length = 0;
total_counter = Counter()
for name, seqrecord in SeqIO.to_dict(SeqIO.parse(args.path, "fasta")).items():
    length, nucl_freq, nucl_counter = get_basic_stat(seqrecord)
    total_length += length
    total_counter.update(nucl_counter);
    print("\t".join([name, str(length)]+nucl_freq))
print("\t".join(["total", str(total_length)] + ["%.3f" % (total_counter.get(x, 0)/float(total_length)) for x in NUCLS]))
    