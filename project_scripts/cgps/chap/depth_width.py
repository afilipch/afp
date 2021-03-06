#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores relation between GC drops and binding peaks'''

import argparse
import os
import sys
import numpy as np;
import pandas as pd;
from pybedtools import BedTool, Interval
from Bio import SeqIO




parser = argparse.ArgumentParser(description='Explores relation between GC drops and binding peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks");
parser.add_argument('--gcdrops', nargs = '?', required=True, type = str, help = "Path to the GC drops");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");

args = parser.parse_args();


motif = ("TAATAAA", "TAATTAA", "TATTAAA", "TATTTAA");
def check_motif(interval, genome, motif):
    seq = str(genome[interval.start: interval.end])
    for m in motif:
        if m in seq:
            return True;
    else:
        return False


regions = BedTool(args.path);
genome = next(SeqIO.parse(args.genome, 'fasta')).seq
gcdrops_intervals = BedTool(args.gcdrops);        
print(len(gcdrops_intervals))

gcdrops_regions = gcdrops_intervals.intersect(regions, c = True)


binding_motif_width = [len(x) for x in gcdrops_regions if x[6] != '0' and check_motif(x, genome, motif)]
binding_motif_depth = [float(x.score) for x in gcdrops_regions if x[6] != '0' and check_motif(x, genome, motif)]

binding_nonmotif_width = [len(x) for x in gcdrops_regions if x[6] != '0' and not check_motif(x, genome, motif)]
binding_nonmotif_depth = [float(x.score) for x in gcdrops_regions if x[6] != '0' and not check_motif(x, genome, motif)]

nonbinding_motif_width = [len(x) for x in gcdrops_regions if x[6] == '0' and check_motif(x, genome, motif)]
nonbinding_motif_depth = [float(x.score) for x in gcdrops_regions if x[6] == '0' and check_motif(x, genome, motif)]

nonbinding_nonmotif_width = [len(x) for x in gcdrops_regions if x[6] == '0' and not check_motif(x, genome, motif)]
nonbinding_nonmotif_depth = [float(x.score) for x in gcdrops_regions if x[6] == '0' and not check_motif(x, genome, motif)]


import matplotlib.pyplot as plt 
fig, ax = plt.subplots(figsize=(16,9))

plt.plot(binding_nonmotif_depth, binding_nonmotif_width, 'g.', label = 'Bound GC drops without motif')
plt.plot(binding_motif_depth, binding_motif_width, 'b.', label = 'Bound GC drops with motif')    
plt.plot(nonbinding_nonmotif_depth, nonbinding_nonmotif_width, color = 'black', marker = '.', linestyle = '', label = 'Unbound GC drops without motif')
plt.plot(nonbinding_motif_depth, nonbinding_motif_width, color = 'red', marker = '.', linestyle = '', label = 'Unbound GC drops with motif')

plt.xlabel('Minimum GC content', fontsize='x-large')
plt.ylabel('Length of GC drop', fontsize='x-large')
plt.tick_params(axis='both', labelsize='x-large', top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.gca().invert_xaxis()
#print(len(binding_width))
plt.legend(fontsize='x-large', frameon=False)
if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format)
else:
    plt.show()
    
    




