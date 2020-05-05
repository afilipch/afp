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
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the annotated transcripts");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to a folder with statistics");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plots Format");
args = parser.parse_args();


def get_closest(tss, transcript_list, strand):
    if(strand == '+'):
        return min([ x.start - tss.start for x in transcript_list], key = lambda y: abs(y))
    else:
        return min([ tss.start - x.stop + 1 for x in transcript_list], key = lambda y: abs(y))





transcript_dict = defaultdict(list);
for x in BedTool(args.transcripts):
    transcript_dict[(x.chrom, x.strand)].append(x)

tss_dict = defaultdict(list);
for x in BedTool(args.path):
    tss_dict[(x.chrom, x.strand)].append(x)
    
distance_dict = defaultdict(list)
for (chrom, strand), tss_list in tss_dict.items():
    #print(len(tss_list))
    transcript_list = transcript_dict[(chrom, strand)]
    for tss in tss_list:
        distance_dict[(chrom, strand)].append(get_closest(tss, transcript_list, strand));


    
#print(distance_dict)

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
        
    #print(xticklabels, bars)



    fig, ax = plt.subplots(figsize=(16, 9))
    ax.bar(brange, bars, 0.5, color='lightblue')
    plt.xticks(brange, xticklabels)
    plt.ylabel("Fraction of TSS")
    plt.savefig(os.path.join(args.outdir, "%s.%s" % (name, args.format)), format = args.format)
    plt.clf()



for key, distances in distance_dict.items():
    draw_tss_distances(distances, "%s(%s)" % key, DR_TR)
