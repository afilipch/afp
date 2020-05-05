#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts annotated peaks gff file into tsv file'''

import argparse
import sys

from pybedtools import BedTool
from afbio.pybedtools_af import read_comments



parser = argparse.ArgumentParser(description='Converts annotated peaks gff file into tsv file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated peaks, gff");
args = parser.parse_args();

labels  = read_comments(args.path)[0][1:].strip().split(",")
    
attrs = ['gene', 'genesymbol', 'annotation', 'function', 'tss', 'atg', 'gtype', 'anti_gene', 'anti_genesymbol', 'anti_tss']
header = ['Chrom', 'Start', 'End', 'Strand'] + labels + attrs
print("\t".join(header))
for interval in BedTool(args.path):
    print("\t".join( [str(x) for x in (interval.chrom, interval.start, interval.end, interval.strand)] + interval.attrs["topcoverage"].split(',') + [interval.attrs[x] for x in attrs]))
