#! /usr/bin/python
'''Sorts bed/gff file with a respect to strandness'''
import argparse
import sys;
import os;
from collections import defaultdict;

from pybedtools import BedTool;



parser = argparse.ArgumentParser(description='Sorts bed/gff file with a respect to strandness');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the bed/gff file to be sorted");
args = parser.parse_args();


def clear_dict(d):
    for v in d.values():
        v[:]=[];
        
    
def flush(chrom_container, start_container):
    for strand, intervals in start_container.items():
        chrom_container[strand].extend(list(sorted(intervals, key = lambda x: x.end)));
    

sortedbed = BedTool(args.path).sort();
if(len(sortedbed) == 0):
    sys.stderr.write("input file is empty\n");
    sys.exit();



chrom_container = defaultdict(list)
start_container = defaultdict(list)

start_container[sortedbed[0].strand].append(sortedbed[0])
chrom = sortedbed[0].chrom;
start = sortedbed[0].start

for interval in sortedbed[1:]:
    if(chrom == interval.chrom):
        if(start == interval.start):
            start_container[interval.strand].append(interval);
        else:
            flush(chrom_container, start_container);
            clear_dict(start_container)
            start = interval.start;
            start_container[interval.strand].append(interval);
        
    else:
        flush(chrom_container, start_container);
        for strand, intervals in chrom_container.items():
            for si in intervals:
                sys.stdout.write(str(si))
                
        chrom = interval.chrom;
        start = interval.start;
        clear_dict(start_container);
        clear_dict(chrom_container);
        start_container[interval.strand].append(interval);
		
else:		
    flush(chrom_container, start_container);
    for strand, intervals in chrom_container.items():
        for si in intervals:
            sys.stdout.write(str(si))
