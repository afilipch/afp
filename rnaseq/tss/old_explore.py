#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores TSS'''
import argparse
import os
import sys
import numpy as np;
import matplotlib.pyplot as plt;
from collections import defaultdict

from pybedtools import BedTool
from Bio import SeqIO

from afbio.sequencetools import coverage2dict

parser = argparse.ArgumentParser(description='Explores TSS');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the TSS in bed format");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the reference genome, fasta format");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the annotated transcripts");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to a folder with statistics");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plots Format");
args = parser.parse_args();


def transcript_distance(transcripts, covdict, strand):
    distances = [];
    tr_starts = defaultdict(list);
    for tr in transcripts:
        if(tr.strand == strand):
            tr_starts[tr.chrom].append(tr.start)
    for k, coverage in covdict.items():
        starts = tr_starts[k];
        for pos, intensity in enumerate(coverage):
            if(intensity):
                distances.append(min([x-pos for x in starts], key = lambda y: abs(y)))
    return distances;
        
        

def closest_distances(covdict):
    distances = [];
    for coverage in covdict.values():
        curintsity = coverage[0];
        curpos = 0;
        for pos, intensity in enumerate(coverage):
            if(intensity):
                if(curintsity):
                    distances.append(pos-curpos);
                curpos = pos;
                curintsity = intensity;
    return distances;

def remove_neighbors(covdict, maxd):
    res = {}
    for k, coverage in covdict.items():
        selected = []
        for pos, intensity in enumerate(coverage[maxd:], start = maxd):
            if(intensity):
                if(intensity >= max(coverage[pos-maxd: pos+maxd+1])):
                    selected.append(pos);
        newcov = [x[1] if (x[0] in selected) else 0 for x in enumerate(coverage)] 
        res[k] = newcov
    return res;


def stat_cov(covdict):
    allcov  = [];
    for v in covdict.values():
        allcov.extend(v)
    ###Output basic coverage statistics
    sys.stderr.write("###Explore TSS\n")
    signal = [x for x in allcov if x]
    median = np.median(signal)
    sys.stderr.write( "Number of TSS:\t%d\nMedian TSS intensity:\t%1.2f\nMean TSS intensity:\t%1.2f\nMax TSS intensity:\t%1.2f\n\n" % (len(signal), median, np.mean(signal), max(allcov)) )
    return median

covdict = coverage2dict(args.path)
transcripts = BedTool(args.transcripts)
median = stat_cov(covdict)

    
threshold = median*10
covdict = dict([ (x[0], [0 if y<=threshold else y for y in x[1]]) for x in covdict.items() ])
#for k in covdict.values():
    #v = [0 if x<=threshold else x for x in v];
    
stat_cov(covdict)



################# Plotting section ###################################
DR_TSS  = [1 ,2, 3, 5, 10, 20, 50, 100, 200, 500, 1000]
DR_TR  = [-500 , -200, -100, -50, -20, -10, -5, 0, 1, 5, 10, 20, 50, 100, 200, 500, 1000]


def draw_tss_distances(distances, name, drange):
    bars = []
    xticklabels = [];
    norma = len(distances)
    for l, u in zip(drange, drange[1:]):
        if(u-l == 1):
            xticklabels.append(str(l));
            bars.append(len([ x for x in distances if x==l])/norma)
        else:
            xticklabels.append("%d-%d" % (l, u))
            bars.append(len([ x for x in distances if l<=x<=u])/norma)
    xticklabels.append(">%d" % drange[-1])
    bars.append(len([ x for x in distances if x>drange[-1]])/norma)
    brange = range(len(xticklabels))
        
    print(xticklabels, bars)



    fig, ax = plt.subplots(figsize=(16, 9))
    ax.bar(brange, bars, 0.5, color='lightblue')
    plt.xticks(brange, xticklabels)
    plt.ylabel("Fraction of TSS")
    plt.savefig(os.path.join(args.outdir, "%s.%s" % (name, args.format)), format = args.format)
    plt.clf()



distances = closest_distances(covdict);
draw_tss_distances(closest_distances(covdict), "tss_distance_raw", DR_TSS) 
draw_tss_distances(transcript_distance(transcripts, covdict, strand='+'), "transcript_distance_raw", DR_TR)
covdict = remove_neighbors(covdict, 5)
draw_tss_distances(closest_distances(covdict), "tss_distance_filtered", DR_TSS) 
draw_tss_distances(transcript_distance(transcripts, covdict, strand='+'), "transcript_distance_filtered", DR_TR)
