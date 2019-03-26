#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts gff file into tsv file'''

import argparse
import sys

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Converts gff file into tsv file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the gff");
parser.add_argument('--header', nargs = '?', const=True, default=False, type = bool, help = "If set, a header is added to the first line of output tsv file");
parser.add_argument('--attrs', nargs = '+', default='', type = str, help = "GFF attributes to be added to the output tsv file. If not set, all attributes are added");
#parser.add_argument('--maxshift', nargs = '?', default=50, type = int, help = "Max allowed shift (in nucleotides) of the peak top position downstream to start of the gene, to be still counted as peak upstream the gene");
args = parser.parse_args();

def get_attrs(interval, selection):
    if(selection):
        return [interval.attrs[x] for x in selection]
    else:
        return [a[1] for a in sorted(interval.attrs.items(), key = lambda x: x[0]) if a[0] != "Name"];
            
            
def get_header(interval, selection):
    const =  ["Name", "Chrom", "Start", "End", "Strand", "Score"]
    if(selection):
        return const + selection
    else:
        return const + [a[0] for a in sorted(interval.attrs.items(), key = lambda x: x[0]) if a[0] != "Name"];
    
intervals = BedTool(args.path);
if(args.header):
    print("\t".join(get_header(intervals[0], args.attrs)))

for interval in intervals:
    print("\t".join( [str(x) for x in (interval.name, interval.chrom, interval.start, interval.end, interval.strand, interval.score)] + get_attrs(interval, args.attrs) ) )