#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Provides basic description of the input fasta file'''

import argparse
import os
import sys


from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt;
from afbio.numerictools import CDF


parser = argparse.ArgumentParser(description='Provides basic description of the input fasta file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the fasta file");
#parser.add_argument('--genes', nargs = '?', required=True, type = str, help = "Path to the genes, sorted gff format");
parser.add_argument('--verbose', nargs = '?', const=True, default = False, type = bool, help = "Set verbose output");
args = parser.parse_args();

NUCLS = 'ACTGN'
def get_basic_stat(seqrecord):
    return len(seqrecord), Counter(str(seqrecord.seq).upper())
        

if(args.verbose):
    stat_list = [];


length_list = []
nucl_counter = Counter()
for seqrecord in SeqIO.parse(args.path, "fasta"):
    length, nucl = get_basic_stat(seqrecord)
    length_list.append(length)
    nucl_counter.update(nucl);
    if(args.verbose):
        stat_list.append((seqrecord.name, length, nucl))
 
nucl_types =  list(sorted(nucl_counter.keys()))
header = ["Name", "Length"] + nucl_types
print("\t".join(header))
total_length = sum(length_list);
print("\t".join(["total", str(total_length)] + ["%.2f" % ((nucl_counter.get(x, 0)/float(total_length))*100) for x in nucl_types]) )


### VERBOSE ###
if(args.verbose):
    for name, length, nucl in stat_list:
        print("\t".join([name, str(length)] + ["%.2f" % ((nucl.get(x, 0)/float(length))*100) for x in nucl_types]) )
        
        
        
###############################################################################################################    
###Plotting
step = 1000
distances = [x if x<step else step for x in length_list]
xvals, yvals = CDF(distances, zerovalue=0)

fontsize, linewidth = 28, 5

fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])
ax.set_xlabel('sequence length', fontsize=fontsize)
ax.set_ylabel('CDF', fontsize=fontsize)
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for axis in ['bottom','left','right']:
    ax.spines[axis].set_linewidth(linewidth)

ax.plot(xvals, yvals, color = 'darkblue', linewidth=linewidth)
#plt.savefig(os.path.join(args.outdir, "utr5_length.%s"  %  args.format) , format = args.format)
plt.show()
    
