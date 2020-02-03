#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Counts number of AT steps per AT-rich stretch'''

import argparse
import os
import sys
from itertools import product;
import numpy as np;
import pandas as pd;
from pybedtools import BedTool, Interval
from Bio import SeqIO




parser = argparse.ArgumentParser(description='Counts number of AT steps per AT-rich stretch');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks");
parser.add_argument('--atstretches', nargs = '?', required=True, type = str, help = "Path to the AT stretches, gff format");
parser.add_argument('--motif', nargs = '?', default = False, const=True, type = str, help = "If set only the AT stretches with the motif will be considered");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");
parser.add_argument('--matrix', nargs = '?', type = str, help = "Path for the output matrix, tsv format");
args = parser.parse_args();




regions = BedTool(args.path);
atstretches = BedTool(args.atstretches);

#count number of at steps and detect motif

if(args.motif):
    motif = [('A', 'T'), ('A', 'T'), ('A', 'T'), ('A', 'T'), ('A', 'T'), ('A',), ('A', 'T'), ('T'), ('A',), ('A',)]
    #motif = [('A', 'T'), ('A', 'T'), ('A', 'T'), ('A', 'T'), ('A', 'T'), ('A', 'T'), ('A', 'T'), ('A', 'T'), ('A', 'T'), ('A', 'T')]
    motifs = ["".join(x) for x in product(*motif)]
    
    def check_motif(seq):
        for m in motifs:
            if(m in seq):
                return True;
        else:
            return False;
else:
    def check_motif(seq):
        return True;
    

#def count_at_steps(sequence):
    #count = 0;
    #curat = 0;
    #for s in sequence:
        #if(s in 'AT'):
            #curat += 1;
        #if(s in 'GC'):
            #if(curat>1):
                #count += 1;
            #curat = 0;
    #else:
        #if(curat>1):
            #count += 1;        
    #return count;        

def count_at_steps(sequence):
    prev = ''
    count = 0;
    for s in sequence:
        if(s == 'A' and prev == 'T'):
            count += 1;
            prev = '';
            
        elif(s == 'T' and prev == 'A'):
            count += 1;
            prev = '';       
            
        else: 
            prev = s;
    return count;  



atstretches = BedTool([Interval(x.chrom, x.start, x.end, x.name, str(count_at_steps(x.attrs['seq'])), x.strand )  for x in atstretches if check_motif(x.attrs['seq']) ])
#for interval in atstretches:
    #interval[5] = str(count_at_steps(interval.attrs['seq'])); 
    
    
#get a division for the at stretches based on their lengthes
lengthes = np.array([len(x) for x in atstretches]);
l_division = [7,9, 11, 14, 17, 20, 30, 40, 50]# [ int(x) for x in sorted(list(set(np.percentile(lengthes, np.linspace(0, 100, 18)))))];
l_division[-1] += 1;
#gccounts = [int(x.score) for x in atstretches]
#gc_division = [int(x) for x in sorted(set(gccounts))]
ac_division = [0,1,2,3,4,5,6,7,11, max([int(x.score) for x in atstretches])+1]
#print(gc_division)
        



atstretches_regions = atstretches.intersect(regions, c = True)
cmatrix = np.zeros([len(l_division)-1, len(ac_division)-1])
nummatrix = np.zeros([len(l_division)-1, len(ac_division)-1])
for i, (l1, l2) in enumerate(zip(l_division, l_division[1:])):
    for j, (a1, a2) in enumerate(zip(ac_division, ac_division[1:])):
        bound = len([x for x in atstretches_regions if (x[-1] != '0' and len(x) >= l1 and len(x) < l2 and int(x.score) >= a1 and int(x.score) < a2)])
        unbound = len([x for x in atstretches_regions if (x[-1] == '0' and len(x) >= l1 and len(x) < l2 and int(x.score) >= a1 and int(x.score) < a2)])
        if(bound or unbound):
            cmatrix[i,j] = bound/(bound+unbound)
            nummatrix[i,j] = bound+unbound

if(args.matrix):
    np.savetxt(args.matrix, cmatrix, delimiter = "\t", fmt = "%1.3f");
#sys.exit()




import matplotlib.pyplot as plt 
    
xnames = ["%d" % x[0] if x[0] == x[1] - 1 else "%d-%d" % (x[0], x[1]-1) for x in zip(ac_division[:-1], ac_division[1:])]
ynames = ["%d-%d" % (x[0], x[1]-1) for x in zip(l_division[:-1], l_division[1:])]
cmap="GnBu"

fig, ax = plt.subplots(figsize=(9,9))
im = ax.imshow(cmatrix, cmap=cmap, vmax=1)
cbar = ax.figure.colorbar(im, ax=ax, cmap=cmap)
cbar.ax.set_ylabel('Fraction of targeted by CgpS', rotation=-90, va="bottom")

ax.set_xticks(np.arange(len(xnames)))
ax.set_yticks(np.arange(len(ynames)))
ax.set_xticklabels(xnames)
ax.set_yticklabels(ynames)
ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

for edge, spine in ax.spines.items():
    spine.set_visible(False)

ax.set_xticks(np.arange(cmatrix.shape[1]+1)-.5, minor=True)
ax.set_yticks(np.arange(cmatrix.shape[0]+1)-.5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
ax.tick_params(which="minor", bottom=False, left=False)
ax.set_xlabel("Number of A/T steps (at least 2nt)");
ax.set_ylabel("Length of AT-rich stretch");

for i in range(nummatrix.shape[0]):
    for j in range(nummatrix.shape[1]):
        text = ax.text(j, i, "%d" % nummatrix[i, j],
                       ha="center", va="center")

if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format)
else:
    plt.show()
    

