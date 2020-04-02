#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores TSS'''
import argparse
import os
import sys
import numpy as np;
import matplotlib.pyplot as plt;
from collections import defaultdict
from matplotlib.patches import Rectangle, Arrow

from pybedtools import BedTool
from Bio import SeqIO

from afbio.sequencetools import coverage2dict

parser = argparse.ArgumentParser(description='Explores TSS');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the TSS in bed format");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the reference genome, fasta format");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the annotated transcripts");
parser.add_argument('--selected', nargs = '?', required=True, type = str, help = "Path to the selected transcripts");
parser.add_argument('--upstream', nargs = '?', default=1000, type = int, help = "Number of basepairs to take upstream");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to a folder with statistics");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plots Format");
args = parser.parse_args();



transcripts = dict([ (x.attrs['cg'][2:], x) for x in BedTool(args.transcripts) ])
tss_dict = defaultdict(list);
for x in BedTool(args.path):
    tss_dict[(x.chrom, x.strand)].append( (x.start, float(x.score)) )



selected = [];
with open(args.selected) as f:
    for l in f:
        a = l.strip()[2:].split("-")
        if(len(a) == 1):
            selected.append(a)
        else:
            temp = [int(x) for x in a]
            selected.append([str(x) for x in range(min(temp), max(temp)+1) ])
        #print(selected[-1])
        #print(l)

selected = [[transcripts.get(y) for y in x] for x in selected]
#for group in selected:
    #for cg in group:
        #print(cg)
    #print()
    #print("*"*120)


### PLOTTING ###
TR_COLORS = 'limegreen', 'purple'
fontsize, linewidth = 20, 5

def draw_transcript_tss(group, local_tss, x_range, strand, name):
    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])
    ax.set_xlabel('genomic position [bp]', fontsize=fontsize)
    ax.set_ylabel('TSS intensity', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)    
    
    xvals = [x[0] for x in local_tss]
    yvals = [x[1] for x in local_tss]
    
    
    width = max(yvals)/20;
    ytr = -width*1.2
    for tr in group:
        rect = Rectangle( (tr.start, ytr), len(tr), width, facecolor = TR_COLORS[0], edgecolor = TR_COLORS[0] )
        ax.add_patch(rect)
    
    if(strand == '+'):
        for tr in group:
            l = len(tr)/2
            arrow = Arrow(tr.start, ytr + width/2, l, 0, width=width*0.75, facecolor = 'black', edgecolor = 'black')
            ax.text(tr.start + l*1.5, ytr + width/2, tr.attrs['cg'], fontsize = fontsize, color = 'black', horizontalalignment='center', verticalalignment='center')
            ax.add_patch(arrow)
    else:
        for tr in group:
            l = len(tr)/2
            arrow = Arrow(tr.stop, ytr + width/2, -l, 0, width=width*0.75, facecolor = 'black', edgecolor = 'black')
            ax.text(tr.stop - l*1.5, ytr + width/2, tr.attrs['cg'], fontsize = fontsize, color = 'black', horizontalalignment='center', verticalalignment='center')
            ax.add_patch(arrow)
            
    ax.bar(xvals, yvals, color = TR_COLORS[1], width = len(x_range)/400)

    #plt.show()
    
    plt.savefig(os.path.join(args.outdir, "%s.%s"  %  (name, args.format) ) , format = args.format)
    plt.close()
    plt.clf()


print("\t".join(["gene block", "chrom", "start", "stop", "strand", "tss position:intensity"]))   
for group in selected:
    if(all(group)):
        strand = group[0].strand
        chrom = group[0].chrom
        if(strand == '+'):
            x_range = range(min([x.start for x in group]) - args.upstream, max([x.stop for x in group]) + 10)
        else:
            x_range = range(min([x.start for x in group]) - 10, max([x.stop for x in group]) + args.upstream)
        
        local_tss = [x for x in tss_dict[(chrom, strand)] if x[0] in x_range]
        temp = [int(x.attrs['cg'][2:]) for x in group]
        if(len(temp)>1):
            name = "cg%d_%d" % (min(temp), max(temp))
        else: 
            name = group[0].attrs['cg']
                
        draw_transcript_tss(group, local_tss, x_range, strand, name)
        
        
        start = min([x.start for x in group])
        stop = min([x.stop for x in group])
        tss_string= "\t".join(["%d:%d" % x for x in sorted(local_tss, key = lambda x: x[1], reverse = True)])
        print("%s\t%s\t%d\t%d\t%s\t%s" % (name, chrom, start, stop, strand, tss_string))
        
        #with open(os.path.join(args.outdir, 'plots', "%s.tsv"  %  name), 'w') as f:
            #local_tss.sort(key = lambda x: x[1]);
            #header = ["cg", "start", "stop", "strand"] + ["tss_%d" % x[0] for x in local_tss]
        
    
