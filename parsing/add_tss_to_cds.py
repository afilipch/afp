#! /usr/bin/python
'''Adds TSS to the CDS in order to get almost complete genes (without 3utrs)'''

import argparse
import sys
import os
from collections import defaultdict, Counter

import numpy as np;
from pybedtools import BedTool
import matplotlib.pyplot as plt;

from afbio.numerictools import CDF


parser = argparse.ArgumentParser(description='Annotates the provided genomic regions with the distances to the closest transcription start sites');
parser.add_argument('--tss', nargs = '?', required=True, type = str, help = "Path to the TSS annotation file, custom format");
parser.add_argument('--cds', nargs = '?', required=True, type = str, help = "Path to the CDS file, gff format");
parser.add_argument('--distance', nargs = '?', default=300, type = int, help = "Max allowed distance from TSS to the closest CDS");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();

cds_list = BedTool(args.cds)
tss_list = BedTool([x for x in BedTool(args.tss) if x[2] == 'TSS'])
OFFSET = 9;

def parse_closest(interval, offset):
    annotation = interval[offset+8]
    score = float(interval[8].split("_")[-1])
    if(annotation == '.'):
        return [None, interval.start, interval.strand, score, 0]
    else:
        name = dict([x.split("=") for x in annotation.split(";")])['ID']
        distance = int(interval[-1]);
        return [name, interval.start, interval.strand, score, distance]
    
    

annotated_tss = [];
for interval in tss_list.closest(b=cds_list, s = True, iu = True, D='a', io=True):
    annotated_tss.append(parse_closest(interval, OFFSET))

orphans = [x for x in annotated_tss if not x[0]]
orphans_distance = [x for x in annotated_tss if x[0] and x[4] > args.distance]
valid_tss = [x for x in annotated_tss if x[0] and x[4] <= args.distance]

geneid2tss = defaultdict(list);
for vtss in valid_tss:
    geneid2tss[vtss[0]].append(vtss[1]);
    
for cds in cds_list:
    tss_vars = geneid2tss.get(cds.name)
    if(tss_vars):
        #sys.stdout.write(str(cds))
        cds.attrs['tss_variants'] = str(len(tss_vars)) 
        for tss in tss_vars:
            if(cds.strand == '+'):
                cds.start = tss;
            else:
                cds.stop = tss;
            sys.stdout.write(str(cds))
        #print()
        #print("_"*180)
        #print()
    else:
        cds.attrs['tss_variants'] = '1'
        sys.stdout.write(str(cds))
        
    
    
sys.stderr.write("\nOrphans tss: %d\nOrphans distance tss: %d\nValid tss: %d\n"  % tuple([len(x) for x in [orphans, orphans_distance, valid_tss]  ]) )
sys.stderr.write("\nnum_tss_per_cds\tnum_of_cds\n")
tss_per_gene = list(sorted(Counter([len(x) for x in geneid2tss.values()]).items(), key = lambda x: x[0]))
for kv in tss_per_gene:
    sys.stderr.write("%d\t%d\n" % kv)
    
    
    
###############################################################################################################    
###Plotting
step = 1000
distances = [x[4] for x in valid_tss + orphans_distance]
distances = [x if x<step else step for x in distances]
xvals, yvals = CDF(distances, zerovalue=0)

fontsize, linewidth = 28, 5

fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])
ax.set_xlabel('5\'UTR length', fontsize=fontsize)
ax.set_ylabel('CDF', fontsize=fontsize)
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for axis in ['bottom','left','right']:
    ax.spines[axis].set_linewidth(linewidth)

ax.plot(xvals, yvals, color = 'darkblue', linewidth=linewidth)
plt.savefig(os.path.join(args.outdir, "utr5_length.%s"  %  args.format) , format = args.format)

    
    
    
    
    
    
    
    
    
    



