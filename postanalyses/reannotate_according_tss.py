#! /usr/bin/python
'''Re-annotates ustranded genomic regions according to the relationship to TSS'''

import argparse
import sys
import os
from collections import defaultdict
#from bisect import bisect_right, bisect_left

import numpy as np;
from pybedtools import BedTool
from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Annotates the provided genomic regions with the distances to the closest transcription start sites');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genomic regions, gff format");
parser.add_argument('--cds', nargs = '?', required=True, type = str, help = "Path to the cds regions, gff format");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts regions, gff format");
parser.add_argument('--inside', nargs = '?', default=200, type = int, help = "Maximum allowed distance to TSS while inside a gene");
parser.add_argument('--maxd', nargs = '?', default=800, type = int, help = "Maximum allowed distance to TSS");
parser.add_argument('--annotation', nargs = '?', required = True, type = str, help = "Path to the annotation file");
#parser.add_argument('--annotation', nargs = '?', required = True, type = str, help = "If set, peaks are trea");
args = parser.parse_args();

def find_closest(peak, starts_plus, stops_minus, cds_dict, ann_dict, d_threshold, inside, if_bed):
    center = int(peak.name)
    distances = [('+', x[0], x[1]-center) for x in starts_plus]
    distances.extend([('-', x[0], center-x[1]+1) for x in stops_minus]);
    distances = [x for x in distances if x[2]>-1*inside]
    strand, gene_name, mindistance = min(distances, key = lambda x: abs(x[2]))
    #print(gene_name)
    if(abs(mindistance) <= d_threshold):
        #gene_name = "_".join(gene_name.split("_")[:-1])
        cds = cds_dict[gene_name]
        if(strand == '+'):
            atg = cds.start - center
        if(strand == '-'):
            atg = center - cds.end + 1
    
        gene_name = gene_name.split("-")[1]
        ann = ann_dict.get(gene_name);
        if(ann):
            genesymbol = ann[0]
            annotation = ann[1]
            function = ann[2]
        else:
            genesymbol = gene_name
            annotation = 'Unknown'
            function = "Unknown"
        
        if(if_bed):
            return construct_gff_interval(peak.chrom, peak.start, peak.stop, 'consensus_region', score='0', strand=strand, source='un', frame='.', attrs=[("Name", peak.name), ("annotation", annotation), ("function", function), ("gene", gene_name), ("genesymbol", genesymbol), ("tss", mindistance), ("atg", abs(atg))])
        else:
            return construct_gff_interval(peak.chrom, peak.start, peak.stop, 'consensus_region', score='0', strand=strand, source='un', frame='.', attrs=[("Name", peak.attrs['Name']), ("maxcov", peak.attrs['maxcov']), ("zscores", peak.attrs['zscores']), ("peakpos", peak.attrs['peakpos']), ("annotation", annotation), ("function", function), ("gene", gene_name), ("genesymbol", genesymbol), ("tss", mindistance), ("atg", abs(atg))])
    else:
        return None

    


ann_dict = {}
with open(args.annotation) as f:
    next(f)
    for l in f:
        a = l.strip().split(";")
        if(len(a) > 9 and a[1] and a[2]):
            ann_dict[a[1]] = a[2], a[8], a[9]


cds_dict = dict([ (x.name, x) for x in BedTool(args.cds) ])
tr_list = BedTool(args.transcripts)
starts_plus = [ (x.name, x.start) for x in tr_list if x.strand == '+']
stops_minus = [ (x.name, x.stop) for x in tr_list if x.strand == '-']



peaks = BedTool(args.path);
if_bed = args.path.split(".")[-1] == 'bed' 
#print(if_bed)
filtered_out = 0;
for peak in peaks:
    newpeak = find_closest(peak, starts_plus, stops_minus, cds_dict, ann_dict, args.maxd, args.inside, if_bed)
    if(newpeak):
        sys.stdout.write(str(newpeak))
        pass;
    else:
        filtered_out += 1;
        
sys.stderr.write("\nTotal peaks: %d\nPeaks passed distance threshold: %d\n"  % (len(peaks), len(peaks) - filtered_out) )
    






    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
