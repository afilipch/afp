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
parser.add_argument('--plot', nargs = '?', required=False, type = str, help = "Path to the output plot");
parser.add_argument('--proteins', nargs = '+', default = [], type = str, help = "If set, only the listed binding proteins are analysed, by default all binding ptoreins are analysed");
parser.add_argument('--numpeaks', nargs = '?', default = False, const = True, type = str, help = "If set, number of peaks inside the region is reported, fraction of binding is reported otherwise");
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



def get_enrichment_fraction(peaklist, region):
    fractions = [];
    for peakfile in peaklist:
        peaks = BedTool(peakfile)
        outside = sum([float(x.attrs['topcoverage']) for x in peaks.intersect(region, f = 0.5, v = True)])
        inside = sum([float(x.attrs['topcoverage']) for x in peaks.intersect(region, f = 0.5, u = True)])    
        if(outside + inside):
            fractions.append(inside/(inside + outside))
    if(fractions):
        return np.mean(fractions)*100, len(fractions)
    else:
        return 'none'
    
    
def get_enrichment_numpeaks(peaklist, region):
    fractions = [];
    for peakfile in peaklist:
        peaks = BedTool(peakfile)
        outside = sum([float(x.attrs['topcoverage']) for x in peaks.intersect(region, f = 0.5, v = True)])
        inside = sum([float(x.attrs['topcoverage']) for x in peaks.intersect(region, f = 0.5, u = True)])    
        if(outside + inside):
            fractions.append(inside)
    if(fractions):
        return int(np.mean(fractions)), len(fractions)
    else:
        return 'none'
    
if(args.numpeaks):
    get_enrichment = get_enrichment_numpeaks
else:
    get_enrichment = get_enrichment_fraction



if(args.proteins):
    labelweigths = [ (x, get_enrichment(protein2exp[x], region)) for x in args.proteins]    
else:
    labelweigths = [ (x[0], get_enrichment(x[1], region)) for x in protein2exp.items() ]
            

 #plot section            
labels = [x[0] for x in labelweigths if x[1] != 'none']
weights = [x[1][0] for x in labelweigths if x[1] != 'none']    
numbers = [str(x[1][1]) for x in labelweigths if x[1] != 'none'] 
print(weights)
            
#prepare an appropriate layout
plt.figure(1, figsize=(16,10), frameon=False)
ax = plt.subplot(111)
plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.96])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
if(args.numpeaks):
    ylim = max(weights)*1.2
else:
    ylim = 100
plt.axis((0, len(labels)+1, 0, ylim))

#set bins and boundaries
boundaries = range(0, len(labels));
bins = range(0, len(labels)+1);

#plot real and control binding pattern
plt.hist(boundaries, weights=weights, bins=bins, align='right', rwidth=0.7, color='skyblue')

if(not args.numpeaks):
    expected = 6.60787254249334
    plt.axhline(y=expected, color='darkblue', linestyle='-')
    plt.yticks([0, expected, 20, 40, 60, 80, 100]);
    ax.set_yticklabels(['0', 'expected\nfraction', '20', '40', '60', '80', '100'])

#set labels and title
if(args.numpeaks):
    plt.ylabel('Total binding [in averaged units]', fontsize = 'x-large')
else:
    plt.ylabel('Fraction of prophage binding [%]', fontsize = 'x-large')
plt.xlabel('Binding protein', fontsize = 'x-large')

#set text labels
for boundary, weight, text in zip(boundaries, weights, numbers):
    ax.text(boundary + 0.85, weight + 5, text, fontsize='x-large');

#set xlabels
plt.xticks(range(1, len(labels)+1));
ax.set_xticklabels(labels, rotation=90)
ax.tick_params(axis='both', labelsize='x-large')

#output plot
if(args.plot):
    plt.savefig(args.plot, bbox_inches='tight', format=args.plot.split(".")[-1])
else: 
	plt.show();           
            
            
            
            
            
            
