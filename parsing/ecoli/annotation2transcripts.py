'''Adds TSS to the CDS extracted from NCBI annotation order to get almost complete transcripts (without 3utrs)'''

import argparse
import sys
import os
from collections import defaultdict, Counter

import numpy as np;
from pybedtools import BedTool
import matplotlib.pyplot as plt;

from afbio.numerictools import CDF


parser = argparse.ArgumentParser(description='Adds TSS to the CDS extracted from NCBI annotation order to get almost complete transcripts (without 3utrs). TSS are considered to be stranded');
parser.add_argument('--tss', nargs = '?', required=True, type = str, help = "Path to the TSS annotation file, custom format");
parser.add_argument('--ncbi', nargs = '?', required=True, type = str, help = "Path to the ncbi annotation");
parser.add_argument('--maxd', nargs = '?', default=500, type = int, help = "Max allowed distance from TSS to the closest CDS");
parser.add_argument('--inside', nargs = '?', default=1, type = int, help = "Max allowed distance downstream from translation start for a TSS");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();


def transcript_length_filter(tss_vars_, cds):
    if(cds.strand == '+'):
        return [x for x in tss_vars_ if cds.stop-x[0]>50]
    else:
        return [x for x in tss_vars_ if x[0]-cds.start>50]

def add_tss(tss, starts_minus, starts_plus, inside):
    if(tss[1]  == '+'):
        temp = [(x[0], x[1] - tss[0]) for x in starts_plus]
    else:
        temp = [(x[0], tss[0] - x[1]) for x in starts_minus]
    temp = [x for x in temp if x[1] >= -1*inside]
    if(temp):
        return min(temp, key = lambda x: abs(x[1]))

tss_list = []
with open(args.tss) as f:
    next(f)
    for l in f:
        a = l.strip().split("\t")
        if(len(a) == 4):
            stationary = 0;
        else:
            stationary = 1;
        log = int(bool(a[3]))
        tss_list.append(( int(a[0])-1, a[1], log, stationary))
        
cds_list = [x for x in BedTool(args.ncbi) if x[2] == 'CDS']
starts_plus = [((x.name, x.start), x.start) for x in cds_list if x.strand == '+']
starts_minus = [((x.name, x.start), x.stop) for x in cds_list if x.strand == '-']

name2info = defaultdict(list)
distances = [];
orphans = 0;
for tss in tss_list:
    res = add_tss(tss, starts_minus, starts_plus, args.inside)
    if(res):
        name, distance = res
        distances.append(distance);
        if(distance <= args.maxd):
            name2info[name].append(( tss[0], distance, tss[2], tss[3] ))
    else:
        orphans+=1

sys.stderr.write("\nTotal tss: %d\nValid tss: %d\nOrphan tss: %d\n"  % (len(distances), sum([len(x) for x in name2info.values()]), orphans ) )
sys.stderr.write("\nnum_tss_per_cds\tnum_of_cds\n")
tss_per_gene = list(sorted(Counter([len(x) for x in name2info.values()]).items(), key = lambda x: x[0]))
for kv in tss_per_gene:
    sys.stderr.write("%d\t%d\n" % kv)
    

for cds in cds_list:
    cds.attrs['cds'] = "%d:%d" % (cds.start, cds.stop)
    tss_vars = transcript_length_filter(name2info[(cds.name, cds.start)], cds)
    if(tss_vars):
        for tss in tss_vars:
            if(cds.strand == '+'):
                cds.start = tss[0];
            else:
                cds.stop = tss[0];
                
            cds.attrs['tss_variants'] = str(len(tss_vars))
            cds.attrs['tss_detected'] = '1'
            cds.attrs['tss_distance'] = str(tss[1])
            cds.attrs['log'] = str(tss[2])
            cds.attrs['stationary'] = str(tss[3])
            sys.stdout.write(str(cds))
    else:
        cds.attrs['log'] = 'unknown'
        cds.attrs['stationary'] = 'unknown'
        cds.attrs['tss_variants'] = '1'
        cds.attrs['tss_detected'] = '0'
        cds.attrs['tss_distance'] = '0'
        sys.stdout.write(str(cds))




#sys.exit();
###############################################################################################################    
###Plotting
step = 1000
#distances = [x[3] for x in valid_tss + orphans_distance]
distances = [x for x in distances if x>=0]
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
#plt.show()
