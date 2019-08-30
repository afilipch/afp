#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts bed to gff file'''

import argparse
from pybedtools import BedTool



parser = argparse.ArgumentParser(description='Converts bed to gff file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the gff file");
args = parser.parse_args();

for interval in BedTool(args.path):
    print("%s\t%d\t%d\t%s\t%s\t%s" % (interval.chrom, interval.start, interval.stop, interval.name, interval.score, interval.strand));
