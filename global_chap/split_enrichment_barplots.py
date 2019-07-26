#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Checks binding enrichment inside the given locus for multiple binding factors, creates barplots'''

import argparse
import os
import sys
from os import listdir
from os.path import isfile
import numpy as np;
from collections import defaultdict

from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt





parser = argparse.ArgumentParser(description='Checks binding enrichment inside the given locus for multiple binding factors, creates barplots');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with binding peaks. \'peaks\' folder in chipchap project folder");
parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Path to experiments annotation table");
parser.add_argument('--start', nargs = '?', required=True, type = int, help = "Start of the region");
parser.add_argument('--end', nargs = '?', required=True, type = int, help = "End of the region");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory for the plots");
parser.add_argument('--plotformat', nargs = '?', default = 'png', type = str, help = "Plot format");
args = parser.parse_args();

region = BedTool([Interval('chr1', args.start, args.end, strand = '+', score = '0', name = 'region')])

exp2protein = {};
with open(args.table) as f:
    next(f);
    for l in f:
        a = l.strip().split("\t")
        exp2protein[".".join(a[1].split(".")[:-1])] = a[3]
#print(exp2protein)
    
peakfiles = [os.path.join(args.path, f) for f in listdir(args.path) if isfile(os.path.join(args.path, f)) and 'annotated' in f]
protein2exp = defaultdict(list)
for pf in peakfiles:
    name = ".".join( os.path.basename(pf).split(".")[:-2] )
    protein2exp[exp2protein[name]].append(pf)    



def get_enrichment(peaklist, region):
    fractions = [];
    labels = [];
    for peakfile in peaklist:
        peaks = BedTool(peakfile)
        outside = sum([float(x.attrs['topcoverage']) for x in peaks.intersect(region, f = 0.5, v = True)])
        inside = sum([float(x.attrs['topcoverage']) for x in peaks.intersect(region, f = 0.5, u = True)])    
        if(outside + inside):
            fractions.append(inside/(inside + outside))
            labels.append(".".join( os.path.basename(peakfile).split(".")[:-2] ))

    return [x*100 for x in fractions], labels

    
  
  
def barplot(weights, labels, path, format_):
#prepare an appropriate layout
    plt.figure(1, figsize=(16,10), frameon=False)
    ax = plt.subplot(111)
    plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.96])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ylim = 100
    plt.axis((0, len(labels)+1, 0, ylim))

    #set bins and boundaries
    boundaries = range(0, len(labels));
    bins = range(0, len(labels)+1);

    #plot real and control binding pattern
    plt.hist(boundaries, weights=weights, bins=bins, align='right', rwidth=0.7, color='skyblue')

    expected = 6.60787254249334
    plt.axhline(y=expected, color='darkblue', linestyle='-')
    plt.yticks([0, expected, 20, 40, 60, 80, 100]);
    ax.set_yticklabels(['0', 'expected\nfraction', '20', '40', '60', '80', '100'])

    #set labels and title

    plt.ylabel('Fraction of prophage binding [%]', fontsize = 'x-large')
    plt.xlabel('Experiment name', fontsize = 'x-large')


    #set xlabels
    for x, label in enumerate(labels):
        ax.text(x+1, 75, label, rotation = 'vertical', fontsize = 'large')
    #plt.xticks(range(1, len(labels)+1));
    #ax.set_xticklabels(labels, rotation = 90)
    ##plt.rcParams["axes.axisbelow"] = True
    #ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='both', labelsize='large')

    plt.savefig(path, format=format_)
    plt.clf()
          
            
            
for protein, peaklist in protein2exp.items():
    weights, labels = get_enrichment(peaklist, region)
    path = os.path.join(args.outdir, "%s.barplot.%s" % (protein, args.plotformat));
    if(len(weights) >= 5):
        barplot(weights, labels, path, args.plotformat)
            
            
            
