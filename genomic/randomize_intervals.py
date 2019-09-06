import argparse
import os
import sys
import numpy as np;
from collections import defaultdict
from pybedtools import BedTool, Interval
import random




parser = argparse.ArgumentParser(description='Constructs random intervals based on the provided real ones');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks, gff/bed format");
parser.add_argument('--number', nargs = '?', default=1000, type = int, help = "Number of the constructed random intervals");
args = parser.parse_args();

intervals = BedTool(args.path)
chroms = [x.chrom for x in intervals]
chrom2interval = defaultdict(list);
for interval in intervals:
    chrom2interval[interval.chrom].append(interval)
    
    
    
strands, starts, stops, lengths = {}, {}, {}, {} 
for chrom, l_intervals in chrom2interval.items():
    strands[chrom] = [x.strand for x in l_intervals]
    lengths[chrom] = [len(x) for x in l_intervals]
    starts[chrom] = min([x.start for x in l_intervals])
    stops[chrom] = max([x.start for x in l_intervals])
 
#print(chroms, starts, stops)
    
for n in range(1, args.number+1):
    chrom = random.choice(chroms)
    strand = random.choice(strands[chrom])
    length = random.choice(lengths[chrom])
    start = random.randint(starts[chrom], stops[chrom])
    sys.stdout.write(str(Interval(chrom, start, start+length, name = "random_%d" % n, score='0', strand=strand)))
    
    

