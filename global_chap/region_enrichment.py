#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Checks binding enrichment inside the given locus'''

import argparse
import os
import sys
from os import listdir
from os.path import isfile
import numpy as np;

from pybedtools import BedTool, Interval





parser = argparse.ArgumentParser(description='Checks binding enrichment inside the given locus');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with binding peaks. \'peaks\' folder in chipchap project folder");
parser.add_argument('--start', nargs = '?', required=True, type = int, help = "Start of the region");
parser.add_argument('--end', nargs = '?', required=True, type = int, help = "End of the region");
parser.add_argument('--length', nargs = '?', required=True, type = int, help = "Length of the genome");
args = parser.parse_args();

region = BedTool([Interval('chr1', args.start, args.end, strand = '+', score = '0', name = 'region')])
peakfiles = [os.path.join(args.path, f) for f in listdir(args.path) if isfile(os.path.join(args.path, f)) and 'annotated' in f]

def get_enrichment(peakfile, region, length):
    peaks = BedTool(peakfile)
    ilen = len(region[0])
    olen = length - ilen
    outside = [float(x.attrs['topcoverage']) for x in peaks.intersect(region, f = 0.5, v = True)]
    inside = [float(x.attrs['topcoverage']) for x in peaks.intersect(region, f = 0.5, u = True)]
    
    normed_count = (len(inside)/ilen)/(len(outside)/olen)
    normed_sum = (sum(inside)/ilen)/(sum(outside)/olen)
    overrepresented = int(normed_count > 2 or normed_sum > 2)
    return "%s\t%d\t%d\t%1.1f\t%1.1f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%d"  % (os.path.basename(peakfile).split(".")[0] , len(inside), len(outside), sum(inside), sum(outside), np.mean(inside), np.mean(outside), np.median(inside), np.median(outside), normed_count, normed_sum, overrepresented)
    #for interval in peaks.intersect(region, f = 0.5, u = True):
        #sys.stdout.write(str(interval))


#for s in range(0, 4000000, 50000):
    #e = s + 100000
    #region = BedTool([Interval('chr1', s, e, strand = '+', score = '0', name = 'region')])
    #local_coverage, total_coverage = get_enrichment(peakfiles[0], region);
    #print(s, e, local_coverage);
    
print("experiment\tnumber peaks inside\tnumber peaks outside\ttotal peak coverage inside\ttotal peak coverage outside\tmean peak coverage inside\tmean peak coverage outside\tmedian peak coverage inside\tmedian peak coverage outside\tpeak density ration\tpeak total coverage ratio\toverrepresented")
for peakfile in peakfiles:
    print(get_enrichment(peakfile, region, args.length));
