#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Annotates the mappings (bed format) and outputs their distribution among genomic feature types (coding genes, rRNA, tRNA and so on)'''

import argparse
import os
import sys
from collections import defaultdict

import numpy as np;
import matplotlib.pyplot as plt;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Annotates the mappings (bed format) and outputs their distribution among genomic feature types (coding genes, rRNA, tRNA and so on)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the mapping, bed format");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the annotation file, gff format");
parser.add_argument('--multi', nargs = '?', default='', type = str, help = "Path to the multimappers");
parser.add_argument('--stranded', nargs = '?', default=False, const = True, type = str, help = "If set the RNA-seq data are supposed to be stranded");

parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();

def parse_intersection(interval, offset):
    a = interval[offset:];
    if(a[0] == '.'):
        gtype = 'intergenic';
        gname = 'None';
    else:
        #print(a)
        d = dict([tuple(x.strip().split('=')) for x in a[8].split(';')])
        gtype = d['type'];
        gname = d['Name'];
    return gname, gtype
        #print(a)
        #print(gtype, gname)
        #print()
        

def annotate_intervals(intervals, offset, norm, gene_types):
    local_normed = 1.0/(norm*len(intervals));
    for interval in intervals:
        gname, gtype = parse_intersection(interval, offset);
        gene_types[gtype] += local_normed


def getname(interval, multi):
    if(multi):
        a = interval.name.split('_');
        return "_".join(a[:-1]), int(a[-1]);
    else:
        return interval.name, 1



def fill_gene_types(mappings, annotation, gene_types, multi, stranded):
    offset = len(mappings[0].fields)
    curname = '';
    curints = [];
    for interval in mappings.intersect(annotation, s = stranded, wao = True, f = 0.5):
        name, norm = getname(interval, multi);
        if(name == curname):
            curints.append(interval);
        else:
            if(curints):
                #print()
                annotate_intervals(curints, offset, norm, gene_types)
            curints = [interval]
            curname = name
        #print(name, norm);
    else:
        annotate_intervals(curints, offset, norm, gene_types)
        
        
gene_types = defaultdict(float)
mappings = BedTool(args.path)
annotation = BedTool(args.annotation)

fill_gene_types(mappings, annotation, gene_types, False, args.stranded);

if(args.multi):
    mappings = BedTool(args.multi)
    fill_gene_types(mappings, annotation, gene_types, True, args.stranded);
    
    
for kv in sorted(gene_types.items(), key = lambda x: x[0]):
    print("%s:\t%d" % kv)
    
    
### Plot a piechart
raw_labels = list(sorted(gene_types.keys()))
raw_sizes = [gene_types[x] for x in raw_labels]
sizes, labels = [], []
ss, sl = 0, ''
norm = sum(raw_sizes);
threshold = 0.01
for s, l in zip(raw_sizes, raw_labels):
    if(s/norm > threshold):
        sizes.append(s);
        labels.append(l);
    else:
        ss += s;
        sl += '%s (%.2f%%)\n' % (l, 100*s/norm)
if(ss):
    sizes.append(ss)
    labels.append(sl)
        
colordict = {'CDS': 'gold', 'intergenic': 'yellowgreen', 'ncRNA': 'lightcoral', 'rRNA': 'lightskyblue', 'tRNA': 'crimson', 'tmRNA': 'indigo', 'upstream': 'orange',  'other': 'violet'}        
colors = [];
for l in labels:
    colors.append(colordict.get(l, colordict['other']))
        
        
plt.figure(figsize=(12,9))
plt.pie(sizes, explode=None, labels=labels, colors=colors,
    autopct='%1.1f%%', shadow=False, startangle=30, pctdistance=0.8)
plt.axis('equal')


name = os.path.basename(args.path).split(".")[0]
plt.title("experiment: %s\ntotal reads mapped: %d" % (name, norm))
plt.savefig(os.path.join(args.outdir, "%s.annotation.%s" % (name, args.format)), format=args.format)



