#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores relation between AT stretches and binding peaks'''

import argparse
import os
import sys
import numpy as np;
import pandas as pd;
from pybedtools import BedTool, Interval
from Bio import SeqIO




parser = argparse.ArgumentParser(description='Explores relation between AT stretches and binding peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks");
#parser.add_argument('--coverage', nargs = '?', required=True, type = str, help = "Path to the binding coverage track");
parser.add_argument('--atstretches', nargs = '?', required=True, type = str, help = "Path to the AT stretches, bed format");
#parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");

args = parser.parse_args();
#coverage = pd.read_csv(args.coverage, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values

regions = BedTool(args.path);
atstretches = BedTool(args.atstretches);

#get a division for the at stretches based on their lengthes
lengthes = np.array([len(x) for x in atstretches]);
l_division = [7,9, 11, 14, 17, 20, 30, 40, 50]# [ int(x) for x in sorted(list(set(np.percentile(lengthes, np.linspace(0, 100, 18)))))];
l_division[-1] += 1;
gccounts = [int(x.score) for x in atstretches]
gc_division = [int(x) for x in sorted(set(gccounts))]
        



atstretches_regions = atstretches.intersect(regions, c = True)
cmatrix = np.zeros([len(l_division)-1, len(gc_division)])
nummatrix = np.zeros([len(l_division)-1, len(gc_division)])
for i, (l1, l2) in enumerate(zip(l_division, l_division[1:])):
    for j, gc in enumerate(gc_division):
        bound = len([x for x in atstretches_regions if (x[-1] != '0' and len(x) >= l1 and len(x) < l2 and int(x.score) == gc)])
        unbound = len([x for x in atstretches_regions if (x[-1] == '0' and len(x) >= l1 and len(x) < l2 and int(x.score) == gc)])
        if(bound or unbound):
            cmatrix[i,j] = bound/(bound+unbound)
            nummatrix[i,j] = bound+unbound

print(cmatrix);
#sys.exit()




import matplotlib.pyplot as plt 
    
xnames = [str(x) for x in gc_division]
ynames = ["%d-%d" % (x[0], x[1]-1) for x in zip(l_division[:-1], l_division[1:])]
cmap="GnBu"

fig, ax = plt.subplots()
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
ax.set_xlabel("Number of G/C nucleotides");
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
    
    
    
#binding_length = [len(x) for x in atstretches_regions if x[6] != '0']
#binding_gccount = [float(x.score) for x in atstretches_regions if x[6] != '0']
#nonbinding_length = [len(x) for x in atstretches_regions if x[6] == '0']
#nonbinding_gccount = [float(x.score) for x in atstretches_regions if x[6] == '0']

#plt.plot(binding_gccount, binding_length, 'b.', label = 'Bound GC drops')
#plt.plot(nonbinding_gccount, nonbinding_length, 'r.', label = 'Unbound GC drops')
#plt.xlabel('Number of G/C in AT stretch')
#plt.ylabel('Length of AT stretch')
#plt.gca().invert_xaxis()
##print(len(binding_width))
#plt.legend()
#plt.show()



