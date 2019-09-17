'''here we try to find the properties based on AT content which discriminate binding peaks from the other genome'''
import argparse
import os
import sys
import copy
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
parser.add_argument('--clear', nargs = '?', default=100, type = int, help = "Additional flank length around the peak (after adding pflank) to remove from the no-peak regions");
parser.add_argument('--at_flank', nargs = '?', default=1000, type = int, help = "Flank length to calculate AT around the peaks");

parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");

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
regions_dict = dict([ (x.name, x) for x in regions ])

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
        
        
def split_at_track(at_dict, regions, clear):
    hollow_at_dict = copy.deepcopy(at_dict)
    inside = []
    for region in regions:
        score =  max(at_dict[region.chrom][region.start:region.stop])
        inside.append(Interval(region.chrom, region.start, region.end, region.name, str(score), region.strand));
        hollow_at_dict[region.chrom][region.start-clear:region.stop+clear] = 0;
    return inside, hollow_at_dict

def convert_outside(hollow_at_dict, pflank):
    res = [];
    size = pflank*2+1
    for chrom, arr in hollow_at_dict.items():
        for start in list(range(size, len(arr), size*2)):
            end = start+size
            score = max(arr[start:end])
            res.append(Interval(chrom, start, end, str((end+start)//2), str(score), '+'))
            
        
    #print(len(res))
    return res;



def compare_inside_outside(inside, outside):
    percentiles = [0, 10, 20, 30, 40, 50, 60, 70]
    percentiles = [40, 50, 60]
    inside_maxes = [float(x.score) for x in inside]
    #inside_maxes.sort();
    outside_maxes = [float(x.score) for x in outside]
    #outside_maxes.sort();
    #print(inside_maxes)
    ##print(outside_maxes[-60:])
    
    res = []
    for p in percentiles:
        at_threshold = np.percentile(inside_maxes, p) - 0.0001;
        res.append((p, at_threshold*100, len([x for x in inside_maxes if x>=at_threshold]), len([x for x in outside_maxes if x>=at_threshold])*2))
    return res
        
        #at_threshold_1 = np.percentile(outside_maxes, p);
        #print(p, at_threshold, at_threshold_1)
        
        
def get_at_flanks(interval, genome, flank):
    f1 = get_at_content(genome[interval.chrom][interval.start-flank:interval.start-100].seq)
    f2 = get_at_content(genome[interval.chrom][interval.stop+100:interval.stop+flank].seq)
    return (f1+f2)/2

def get_at_flanks_max(interval, genome, flank, window):
    size = window*2+1
    f1 = max([ get_at_content(x) for x in sliding_window(genome[interval.chrom][interval.start-flank:interval.start].seq, size) ])
    f2 = max([ get_at_content(x) for x in sliding_window(genome[interval.chrom][interval.stop:interval.stop+flank].seq, size) ])
    return (f1+f2)/2
        
        
        

        
for window in range(*args.windows):
    at_dict = get_at_dict(genome, window, args.mask)
    inside, hollow_at_dict = split_at_track(at_dict, regions, args.clear);
    outside = convert_outside(hollow_at_dict, args.pflank)
    separation = compare_inside_outside(inside, outside)
    print("#"*180)
    print("AT motif length %d" % (window*2+1))
    for st in separation:
        print("%d\t%1.1f\t%d\t%d" % st)
    print()
    sys.stderr.write("AT motif of length %d has been processed\n" % (window*2+1))
    
    #flank_at_list = []
    #for oreg in [x for x in inside if float(x.score)>0.65]:
        #flank_at = get_at_flanks_max(oreg, genome, args.at_flank, window)
        #flank_at_list.append(flank_at);
        #oreg.chrom = 'chr1'
        #sys.stdout.write("%s\t%1.3f\n" % (str(oreg).strip(), flank_at));
    #print(np.percentile(flank_at_list, 25))
    #print(np.percentile(flank_at_list, 75))
      
    #flank_at_list = [] 
    #for oreg in [x for x in outside if float(x.score)>0.65]:
        #flank_at = get_at_flanks_max(oreg, genome, args.at_flank, window)
        #flank_at_list.append(flank_at);
        #oreg.chrom = 'chr1'
        #sys.stdout.write("%s\t%1.3f\n" % (str(oreg).strip(), flank_at));
    #print(np.percentile(flank_at_list, 25))
    #print(np.percentile(flank_at_list, 75))
    
    
    
