'''here we try to find the properties based on AT content which discriminate binding peaks from the other genome'''
import argparse
import os
import sys
import numpy as np;
from collections import defaultdict, Counter
from pybedtools import BedTool, Interval
from Bio import SeqIO
import matplotlib.pyplot as plt;

from afbio.sequencetools import get_at_content, sliding_window, coverage2dict

#from afbio.sequencetools import get_at_content, sliding_window, array2fixed_length, coverage2dict




parser = argparse.ArgumentParser(description='Checks the AT content along the regions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks, gff format");
parser.add_argument('--windows', nargs = 2, default=(5,20), type = int, help = "Range of the windows sizes (flanks around a particular genomic position) to look for an AT content");
parser.add_argument('--pflank', nargs = '?', default=60, type = int, help = "Length of the flank around a peak center to look for the maximum AT content");
parser.add_argument('--mask', nargs = '?', default=200, type = int, help = "Tails of the chromosomes will be masked with the provided length, mask must be higher than window");
#parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
#parser.add_argument('--zscore', nargs = '?', default=2, type = int, help = "Z-score threshold for the binding peaks");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
args = parser.parse_args();

def convert_region(region, pflank):
    top = (region.end+region.start)//2
    return Interval(region.chrom, top-pflank, top+pflank+1, region.name, region.score, region.strand); 


genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
regions = BedTool(args.path)
regions = [convert_region(x, args.pflank) for x in regions]
#regions = BedTool([x for x in regions if float(x.score)>args.zscore])


def get_at_dict(genome, window, mask):
    at_dict = {};
    masked = mask - window;
    frame = window*2+1
    for chrom, seq in genome.items():
        at = [get_at_content(x) for x in sliding_window(seq[masked:-masked], frame)]
        at = [0]*mask + at + [0]*mask
        at_dict[chrom] = np.array(at);
        #print(at[:20])
        #print(str(seq[:30].seq))
    return at_dict
        
        
def split_at_track(at_dict, regions):
    inside = []
    for region in regions:
        inside.append(list(at_dict[region.chrom][region.start:region.stop]));
        at_dict[region.chrom][region.start:region.stop] = 0;
    return inside, list(at_dict.values())

def convert_outside(outside, pflank):
    res = [];
    for arr in outside:
        size = pflank*2+1
        split = list(range(size, len(arr), size))
        selected = [ x[1] for x in enumerate(np.split(arr, split)) if x[0] % 2]
        res.extend(selected);
    return res;



def compare_inside_outside(inside, outside):
    percentiles = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    inside_maxes = [max(x) for x in inside]
    inside_maxes.sort();
    outside_maxes = [max(x) for x in outside]
    outside_maxes.sort();
    
    for p in percentiles:
        at_threshold = np.percentile(inside_maxes, p);
        at_threshold_1 = np.percentile(outside_maxes, p);
        print(p, at_threshold, at_threshold_1)
        
        
        
        
for window in range(*args.windows):
    at_dict = get_at_dict(genome, window, args.mask)
    inside, outside = split_at_track(at_dict, regions);
    outside = convert_outside(outside, args.pflank)
    compare_inside_outside(inside, outside)
    
