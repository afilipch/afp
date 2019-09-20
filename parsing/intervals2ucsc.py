#! /usr/bin/python
'''Converts provided intervals into ucsc compatible format'''

import argparse
import sys
import os
from collections import defaultdict, Counter

import numpy as np;
from pybedtools import BedTool, Interval

#from afbio.numerictools import CDF


parser = argparse.ArgumentParser(description='Converts provided intervals into ucsc compatible format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the intervals, bed format");
parser.add_argument('--name', nargs = '?', required=True, type = str, help = "Name of the dataset")
parser.add_argument('--description', nargs = '?', required=True, type = str, help = "Description of the dataset")
parser.add_argument('--minscore', nargs = '?', default=0, type = int, help = "Score threshold for the binding peaks");
args = parser.parse_args();

trackopts = "track name=\"%s\" description=\"%s\" visibility=1 color=0,60,120 useScore=1" % (args.name, args.description)
print(trackopts);

regions = BedTool(args.path)
regions = BedTool([x for x in regions if float(x.score)>args.minscore])


score_range = [166, 277, 388, 499, 611, 722, 833, 945, 1000]
percentile_range = list(np.linspace(0,100, len(score_range)+1))
scores = [float(x.score) for x in regions]
thresholds = np.percentile(scores, percentile_range);
thresholds[-1] *= 1.01

updated_regions = [];
for region in regions:
    for ucsc_score, t1, t2 in zip(score_range, thresholds, thresholds[1:]):
        if(t1 <= float(region.score) < t2):
            updated_regions.append(Interval('chr1', region.start, region.end, region.name, str(ucsc_score-2), region.strand))

            
#print([float(x.score) for x in updated_regions])           


for interval in updated_regions:
    sys.stdout.write(str(interval))
