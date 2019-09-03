#! /usr/bin/python
'''Analyses TSS differences on phage versus non-phage'''

import argparse
import sys
import os
from collections import defaultdict, Counter
import matplotlib.pyplot as plt;
import numpy as np;

from pybedtools import BedTool

#from afbio.numerictools import CDF




parser = argparse.ArgumentParser(description='Analyses TSS differences on phage versus non-phage');
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts coordinates augmented with tss information [output of ~/afp/genomic/at_relative_to_transcripts.py], bed/gff file");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();

transcripts = BedTool(args.transcripts);
phages = BedTool(args.phages);
phaged_transcripts = [transcripts.intersect(b=phages, u =True, f = 0.5), transcripts.intersect(b=phages, v =True, f = 0.5)]

tss_counts_list = [];


for ptr in phaged_transcripts:
    temp_dict = defaultdict(list);
    for transcript in ptr:
        temp_dict[transcript.name].append(int(transcript.attrs['tss_variants']))
    tss_counts = [x[0] for x in temp_dict.values()]
    tss_counts = [x if x<5 else 5 for x in tss_counts]
    norma = len(tss_counts)
    tss_counts = Counter(tss_counts);
    tss_counts = [tss_counts[x]/norma for x in range(1, 6)]
    tss_counts_list.append(tss_counts);
    
for x in zip(range(1, 6), tss_counts_list[0], tss_counts_list[1]):
    print("%d\t%1.3f\t%1.3f" % x)

 
############################################################################################################
### Plotting section

### Plot phage versus transcript distribution 

fontsize, linewidth = 28, 5
width = 0.35
labels = [str(x) for x in range(1, 6)]
x = np.arange(len(labels))
barlabels = ['phage', 'non-phage']
barcolors = ['darkblue', 'lightblue']


fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.95])

ax.set_xlabel('Number of TSS per gene', fontsize=fontsize)
ax.set_ylabel('Fraction of genes', fontsize=fontsize)
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for axis in ['bottom','left','right']:
    ax.spines[axis].set_linewidth(linewidth)

ax.bar(x - width/2, tss_counts_list[0], width, label=barlabels[0], color = barcolors[0])
ax.bar(x + width/2, tss_counts_list[1], width, label=barlabels[1], color = barcolors[1])

fig.legend(loc=(0.75, 0.8), frameon=False, fontsize=fontsize, ncol = 1)
plt.savefig(os.path.join(args.outdir, "tss_phage_vs_genomic.%s"  % args.format) , format = args.format)
    
