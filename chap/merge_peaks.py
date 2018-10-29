#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Finds consensus regions for the peaks found in different experiments'''

import argparse
import os
import sys
from collections import defaultdict


import numpy as np;
import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Finds consensus regions for the peaks found in different experiments');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to all the detected peaks");
parser.add_argument('--coverage', nargs = '+', required=True, type = str, help = "Path to the coverage files, must have the order consistent with input \'peak\' files")

parser.add_argument('--zscore', nargs = '?', default=5.0, type = float, help = "Mininimum z-score required for a peak to be a seed for the consensus region");

parser.add_argument('--flank', nargs = '?', default=40, type = int, help = "Flank's length around the detected peaks. Region size = flank*2 + 1");
parser.add_argument('--peakflank', nargs = '?', default=5, type = int, help = "Flank's length around a peak to detect maximum coverage");
args = parser.parse_args();


OVERLAP = 0.55

##########################################################################################################################################
### Read the data and convert them into regions and points 

coverages = [pd.read_csv(x, sep="\t" , names = ["chr", "postion", "coverage"]).coverage.values for x in args.coverage]
peaks_list = [BedTool(x) for x in args.path]

# Processing
        
def peaks2seeds(peaks, flank, zscore_threshold):
    seeds = [(peak.chrom, int(peak.name)-flank, int(peak.name)+flank+1, peak.name, '0', peak.strand) for peak in peaks if float(peak.score)>zscore_threshold];
    return BedTool(seeds)


def peaks2points(peaks):
    return BedTool([(peak.chrom, int(peak.name), int(peak.name)+1, peak.name, peak.score, peak.strand) for peak in peaks])

def getregions(seed_list, overlap):
    regions = seed_list[0]
    for seeds in seed_list[1:]:
        uniq = seeds.intersect(regions, v = True, f = overlap, F = overlap, s = True)
        regions = regions.cat(uniq, postmerge=False);
    return regions;



point_list = [peaks2points(x) for x in peaks_list];
seed_list = [peaks2seeds(x, args.flank, args.zscore) for x in peaks_list];
regions = getregions(seed_list, OVERLAP)



##########################################################################################################################################
### Annotate regions with the information derived from points

def get_best_overlap(intervals):
    if(len(intervals) == 1):
        i = intervals[0]
        if(i[6] == '.'):
            return None;
        else:
            return i[6:]
    else:
        ms = max(float(x[10]) for x in intervals);
        for i in intervals:
            if(float(i[10]) == ms):
                return i[6:]



def get_intersection(points, regions, region_dict):
    curname = '';
    curintervals = [];
    for interval in regions.intersect(points, wao=True, s=True):
        if(interval.name == curname):
            curintervals.append(interval);
        else:
            if(curintervals):
                region_dict[(curintervals[0].chrom, curintervals[0].strand, curname)].append(get_best_overlap(curintervals))
            curintervals = [interval];
            curname =  interval.name
    else:
        region_dict[(curintervals[0].chrom, curintervals[0].strand, curname)].append(get_best_overlap(curintervals))



region_dict = defaultdict(list);
for points in point_list:
    get_intersection(points, regions, region_dict);


##########################################################################################################################################
### Annotate regions with coverage information    
  
def region2gff(region, points, flank, covscores):
    peakpos = ",".join([x[3] if x else 'None' for x in points])
    zscores = ",".join([x[4] if x else 'None' for x in points])
    maxcov = ",".join(covscores);
    return construct_gff_interval(region[0], int(region[2])-flank, int(region[2])+flank+1, 'consensus_region', score='0', strand=region[1], source='un', frame='.', attrs=[('Name', region[2]), ('peakpos', peakpos), ("zscores", zscores), ('maxcov', maxcov)])



def get_maxcov(coverages, points, peakflank):
    covscores = [];
    for an, coverage in zip(points, coverages):
        if(an):
            peakpos = int(an[3])
            covscores.append("%1.3f" % max(coverage[peakpos-peakflank: peakpos+peakflank+1]))
        else:
            covscores.append('None');
    return covscores;

  
for c, (region, points) in enumerate(region_dict.items()):
    covscores = get_maxcov(coverages, points, args.peakflank)
    sys.stdout.write(str(region2gff(region, points, args.flank, covscores)))












































##print(len(regions))
#for s in point_list[0]:
    #sys.stdout.write(str(s))