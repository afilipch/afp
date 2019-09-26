#! /usr/bin/python
'''Re-annotates ustranded genomic regions according to the relationship to TSS'''

import argparse
import sys
import os
from collections import defaultdict
from bisect import bisect_right, bisect_left

import numpy as np;
from pybedtools import BedTool
from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Annotates the provided genomic regions with the distances to the closest transcription start sites');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genomic regions, gff format (output of 'postanalyses/annotate_tss.py')");
parser.add_argument('--cds', nargs = '?', required=True, type = str, help = "Path to the cds regions, gff format");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts regions, gff format");
parser.add_argument('--inside', nargs = '?', default=50, type = int, help = "Maximum allowed distance to TSS while inside a gene");
parser.add_argument('--maxd', nargs = '?', default=700, type = int, help = "Maximum allowed distance to TSS");
parser.add_argument('--annotation', nargs = '?', required = True, type = str, help = "If set, the intervals will be renamed and annotated");
args = parser.parse_args();

def find_closest(peak, starts_plus, stops_minus, cds_dict, ann_dict, d_threshold, inside):
    center = int(peak.name)
    distances = [('+', x[0], x[1]-center) for x in starts_plus]
    distances.extend([('-', x[0], center-x[1]+1) for x in stops_minus]);
    distances = [x for x in distances if x[2]>-1*inside]
    strand, gene_name, mindistance = min(distances, key = lambda x: abs(x[2]))
    if(abs(mindistance) <= d_threshold):
        gene_name = "_".join(gene_name.split("_")[:-1])
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
        
        
        return construct_gff_interval(peak.chrom, peak.start, peak.stop, 'consensus_region', score='0', strand=strand, source='un', frame='.', attrs=[("Name", peak.attrs['Name']), ("maxcov", peak.attrs['maxcov']), ("zscores", peak.attrs['zscores']), ("peakpos", peak.attrs['peakpos']), ("annotation", annotation), ("function", function), ("gene", gene_name), ("genesymbol", genesymbol), ("tss", mindistance), ("atg", abs(atg))])
    else:
        #print(peak, strand, gene_name, mindistance)
        #print("#"*150)
        return None

    


ann_dict = {}
with open(args.annotation) as f:
    next(f)
    for l in f:
        a = l.strip().split(";")
        if(len(a) > 9 and a[1] and a[2]):
            ann_dict[a[1]] = a[2], a[8], a[9]
            
            
#print(ann_dict)

#sys.exit()

cds_dict = dict([ (x.name, x) for x in BedTool(args.cds) ])
tr_list = BedTool(args.transcripts)
starts_plus = [ (x.name, x.start) for x in tr_list if x.strand == '+']
stops_minus = [ (x.name, x.stop) for x in tr_list if x.strand == '-']



peaks = BedTool(args.path);
filtered_out = 0;
for peak in peaks:
    newpeak = find_closest(peak, starts_plus, stops_minus, cds_dict, ann_dict, args.maxd, args.inside)
    if(newpeak):
        sys.stdout.write(str(newpeak))
        pass;
    else:
        filtered_out += 1;
        
sys.stderr.write("\nTotal peaks: %d\nPeaks passed distance threshold: %d\n"  % (len(peaks), len(peaks) - filtered_out) )
    







#name2strand = {"tss_downstream_minus_distance": "-", "tss_downstream_plus_distance": "+", "tss_upstream_minus_distance": "-", "tss_upstream_plus_distance": "+"}

#def select_closest(interval, cds_dict, inside, maxd):
    #outside_gene = "tss_downstream_minus_distance", "tss_downstream_plus_distance"
    #inside_gene ="tss_upstream_minus_distance", "tss_upstream_plus_distance";
    #distances = [(x,int(interval.attrs[x])) for x in outside_gene];
    
    #for name in inside_gene:
        #distance = int(interval.attrs[name])
        #if(distance<=inside):
            #distances.append((name, distance))
            
    #min_name, mindistance = min(distances, key=lambda x: x[1])
    #gene = interval.attrs[min_name.replace("distance", "gene")]
    #strand = name2strand[min_name]
    #cds = cds_dict.get(gene)
    #if(cds and mindistance<=maxd):
        #if(strand == '+'):
            #atg = cds.start - (interval.start + interval.stop)//2
        #else:
            #atg = (interval.start + interval.stop)//2 - cds.stop

     
        #res = construct_gff_interval(interval.chrom, interval.start, interval.stop, 'consensus_region', score='0', strand=strand, source='un', frame='.', attrs=[("Name", interval.attrs['Name']), ("maxcov", interval.attrs['maxcov']), ("zscores", interval.attrs['zscores']), ("peakpos", interval.attrs['peakpos']), ("gene", gene), ("tss", mindistance), ("atg", abs(atg)), ("annotation", mindistance) ])
        #return res
    #else:
        #return None
        
        
    
#cds_dict = dict([ (x.name.split("-")[1], x) for x in BedTool(args.cds) ])
#if(args.annotation):
    #ann_dict = dict([ (x.attrs['geneid'].replace("cg", "NCgl"), x) for x in BedTool(args.annotation) ])
    
    


    
#for interval in BedTool(args.path):
    #res = select_closest(interval, cds_dict, args.inside, args.maxd)
    #if(res):
        #if(args.annotation):
            ##print(res.attrs['gene'])
            #ann = ann_dict.get(res.attrs['gene']);
            #if(ann):
                #genesymbol = ann.name;
                #annotation = ann.attrs['annotation']
                #function = ann.attrs['function']
            #else:
                #genesymbol = res.attrs['gene'].replace("NCgl", "cg")
                #annotation = 'Unknown'
                #function = "Unknown"
            #res.attrs['genesymbol'] = genesymbol
            #res.attrs['annotation'] = annotation
            #res.attrs['function'] = function
            
                
        #sys.stdout.write(str(res))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
