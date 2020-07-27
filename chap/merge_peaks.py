#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Finds consensus regions for the peaks found in different experiments/replicates'''


import argparse
import os
import sys
from collections import defaultdict

from pybedtools import BedTool

from afbio.pybedtools_af import construct_gff_interval
from afbio.peaks import find_shared_peaks, shared_peaks_stat_to_string
from afbio.generators import get_only_files


parser = argparse.ArgumentParser(description='Finds consensus regions for the peaks found in different experiments/replicates');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the directory with all the detected peaks");
parser.add_argument('--maxd', nargs = '?', default=60, type = int, help = "Maximal distance allowed between top positions of the peaks");
parser.add_argument('--replicates', nargs = '+', default = [], type = str, help = "If set, provided peaks are considered as having replicates");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
parser.add_argument('--name', nargs = '?', default='all', type = str, help = "Name of the output file");
args = parser.parse_args();


def gather_files(path, replicates, name):
    res = defaultdict(list)
    files = [x for x in get_only_files(path) if "annotated" in x]
    
    for f in files:
        for r in replicates:
            if(os.path.basename(f).startswith(r)):
                res[r].append(f);
    res[name] = files
    return res;
    
    

def print_compiled(compiled, size):
    temp_d = dict(compiled)
    compiled_processed = [temp_d.get(x, None) for x in range(size)]
    
    area_coverage = ",".join([x.attrs['area_coverage'] if x else '0' for x in compiled_processed])
    topcoverage=",".join([x.attrs['topcoverage'] if x else '0' for x in compiled_processed])
    compiled_processed = [x for x in compiled_processed if x]
    
    compiled = [x[1] for x in sorted(compiled, key = lambda x: x[0])]
    pos = int(sum([int(x.name)*float(x.attrs['topcoverage']) for x in compiled_processed])/sum([float(x.attrs['topcoverage']) for x in compiled_processed]))
    start = min([x.start for x in compiled_processed])
    stop = min([x.stop for x in compiled_processed])
    
    consensus = construct_gff_interval(compiled[0].chrom, start, stop, 'consensus', score='0', strand='.', source='.', frame='.', attrs=[('Name', pos), ('topcoverage', topcoverage), ('area_coverage', area_coverage)])
    
    return str(consensus)
    
    #print(consensus)
    #for cs in compiled_processed:
        #print(cs)
    #print()
    #print("*"*150)
    

            
            
def run_accross_chromosome(bedtools_chr, maxd, size):
    stat_counts = []
    marked = [];
    for c, intervals in enumerate(bedtools_chr):
        for interval in intervals:
            marked.append((c, interval))
    marked.sort(key = lambda x: int(x[1].name));
    #for m in marked:
        #print(m[0], m[1].name)
    selections = [];
    current_selection = [marked[0]];
    for m in marked[1:]:
        if(current_selection and int(m[1].name) - int(current_selection[0][1].name) <= maxd):
            current_nums = [x[0] for x in current_selection]
            if(m[0] not in current_nums):
                current_selection.append(m)
        else:
            print_compiled(current_selection, size);
            stat_counts.append(len(current_selection))
            current_selection = [m]
    else:
        print_compiled(current_selection, size)
        stat_counts.append(len(current_selection))
        
    return stat_counts
            

    
    
        
########################################################################################################    
### Execution Section


file_dict = gather_files(args.path, args.replicates, args.name)
for dname, files in file_dict.items():
    blist = [BedTool(x) for x in files]
    size = len(blist)
    res_total, stat_total_counts = find_shared_peaks(blist, args.maxd)

    with open(os.path.join(args.outdir, "%s.gff" % dname), 'w') as f:
        
        f.write("# %s\n" % ",".join([os.path.basename(x).split(".")[0] for x in files]))
        for compiled in res_total:
            f.write(print_compiled(compiled, size))
    
    sys.stderr.write("\n%s\n" % dname);
    sys.stderr.write(shared_peaks_stat_to_string(stat_total_counts, size))

    
























































#import argparse
#import os
#import sys
#from collections import defaultdict


#import numpy as np;
#import pandas as pd;
#from pybedtools import BedTool, Interval
#import matplotlib.pyplot as plt;

#from afbio.pybedtools_af import construct_gff_interval


#parser = argparse.ArgumentParser(description='Finds consensus regions for the peaks found in different experiments');
#parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to all the detected peaks");
#parser.add_argument('--coverage', nargs = '+', required=True, type = str, help = "Path to the coverage files, must have the order consistent with input \'peak\' files")

#parser.add_argument('--zscore', nargs = '?', default=5.0, type = float, help = "Mininimum z-score required for a peak to be a seed for the consensus region");

#parser.add_argument('--flank', nargs = '?', default=40, type = int, help = "Flank's length around the detected peaks. Region size = flank*2 + 1");
#parser.add_argument('--peakflank', nargs = '?', default=5, type = int, help = "Flank's length around a peak to detect maximum coverage");
#args = parser.parse_args();


#OVERLAP = 0.55

###########################################################################################################################################
### Read the data and convert them into regions and points 

#coverages = [pd.read_csv(x, sep="\t" , names = ["chr", "postion", "coverage"]).coverage.values for x in args.coverage]
#peaks_list = [BedTool(x) for x in args.path]

## Processing
        
#def peaks2seeds(peaks, flank, zscore_threshold):
    #seeds = [(peak.chrom, max(0, int(peak.name)-flank), int(peak.name)+flank+1, peak.name, '0', peak.strand) for peak in peaks if float(peak.score)>zscore_threshold];
    #return BedTool(seeds)


#def peaks2points(peaks):
    #return BedTool([(peak.chrom, int(peak.name), int(peak.name)+1, peak.name, peak.score, peak.strand) for peak in peaks])

#def getregions(seed_list, overlap):
    #regions = seed_list[0]
    #for seeds in seed_list[1:]:
        #uniq = seeds.intersect(regions, v = True, f = overlap, F = overlap, s = True)
        #regions = regions.cat(uniq, postmerge=False);
    #return regions;



#point_list = [peaks2points(x) for x in peaks_list];
#seed_list = [peaks2seeds(x, args.flank, args.zscore) for x in peaks_list];
#regions = getregions(seed_list, OVERLAP)



###########################################################################################################################################
### Annotate regions with the information derived from points

#def get_best_overlap(intervals):
    #if(len(intervals) == 1):
        #i = intervals[0]
        #if(i[6] == '.'):
            #return None;
        #else:
            #return i[6:]
    #else:
        #ms = max(float(x[10]) for x in intervals);
        #for i in intervals:
            #if(float(i[10]) == ms):
                #return i[6:]



#def get_intersection(points, regions, region_dict):
    #curname = '';
    #curintervals = [];
    #for interval in regions.intersect(points, wao=True, s=True):
        #if(interval.name == curname):
            #curintervals.append(interval);
        #else:
            #if(curintervals):
                #region_dict[(curintervals[0].chrom, curintervals[0].strand, curname)].append(get_best_overlap(curintervals))
            #curintervals = [interval];
            #curname =  interval.name
    #else:
        #region_dict[(curintervals[0].chrom, curintervals[0].strand, curname)].append(get_best_overlap(curintervals))



#region_dict = defaultdict(list);
#for points in point_list:
    #get_intersection(points, regions, region_dict);


###########################################################################################################################################
### Annotate regions with coverage information    
  
#def region2gff(region, points, flank, covscores):
    #peakpos = ",".join([x[3] if x else 'None' for x in points])
    #zscores = ",".join([x[4] if x else 'None' for x in points])
    #maxcov = ",".join(covscores);
    #return construct_gff_interval(region[0], max(0, int(region[2])-flank), int(region[2])+flank+1, 'consensus_region', score='0', strand=region[1], source='un', frame='.', attrs=[('Name', region[2]), ('peakpos', peakpos), ("zscores", zscores), ('maxcov', maxcov)])



#def get_maxcov(coverages, points, peakflank):
    #covscores = [];
    #for an, coverage in zip(points, coverages):
        #if(an):
            #peakpos = int(an[3])
            #covscores.append("%1.3f" % max(coverage[max(0, peakpos-peakflank): peakpos+peakflank+1]))
        #else:
            #covscores.append('None');
    #return covscores;

  
#for c, (region, points) in enumerate(region_dict.items()):
    #covscores = get_maxcov(coverages, points, args.peakflank)
    #sys.stdout.write(str(region2gff(region, points, args.flank, covscores)))    
