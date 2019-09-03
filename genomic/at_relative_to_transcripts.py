#! /usr/bin/python
'''Analyses AT content relative to the transcripts and phages positions on a genome'''

import argparse
import sys
import os
from collections import defaultdict
import matplotlib.pyplot as plt;
import numpy as np;

from Bio import SeqIO;
from pybedtools import BedTool;

from afbio.sequencetools import transcript_upstream, upstream_content, transcript_content
from afbio.numerictools import CDF




parser = argparse.ArgumentParser(description='Analyses AT content relative to the transcripts and phages positions on a genome');
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts coordinates, bed/gff file");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--length', nargs = '?', default=40, type = int, help = "Local length of the segment upstream a transcript to check for an AT content");
parser.add_argument('--lookup', nargs = '?', default=500, type = int, help = "Total length of the segment upstream a transcript to check for an AT content");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");

#parser.add_argument('--tss', nargs = '?', required=True, type = str, help = "Path to the TSS annotation file");
#parser.add_argument('--stranded', nargs = '?', const=True, default=False, type = bool, help = "If set, the data considered to be stranded");
args = parser.parse_args();

genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
transcripts = BedTool(args.transcripts);
phages = BedTool(args.phages);


phaged_transcripts = [transcripts.intersect(b=phages, u =True, f = 0.5), transcripts.intersect(b=phages, v =True, f = 0.5)]
transcript_cdf = []
at_content_list = []

for ptr in phaged_transcripts:
    temp_dict = defaultdict(list);
    for transcript in ptr:
        temp_dict[transcript.name].append(transcript_content(transcript, genome))
    at_content = [np.mean(x) for x in temp_dict.values()]
    transcript_cdf.append(CDF(at_content, zerovalue=0))
    at_content_list.append(np.array(at_content))
    
    


upstream_cdf = [];
for ptr in phaged_transcripts:
    at_content = [transcript_upstream(x, genome, 20, 60) for x in ptr];
    at_content = [max(x) for x in at_content if x]
    upstream_cdf.append(CDF(at_content, zerovalue=0))
    at_content_list.append(np.array(at_content))
    
print("type of regions", "num of regions", "mean AT[%]", "median AT[%]", "std AT[%]")
for at_content, name in zip(at_content_list, ["phage transcripts", "non-phage transcripts", "phage upstream", "non-phage upstream"]):
    print("%s\t%d\t%1.1f\t%1.1f\t%1.1f" % (name, len(at_content), at_content.mean()*100, np.median(at_content)*100, at_content.std()*100 ))
    
#print(upstream_cdf[1])
#print(upstream_cdf[0])


############################################################################################################
### Plotting section

### Plot phage versus transcript distribution

def draw_cdf(cdf_list, name, fontsize, linewidth):
    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

    ax.set_xlabel('AT-content', fontsize=fontsize)
    ax.set_ylabel('CDF', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)

    ax.plot(cdf_list[0][0], cdf_list[0][1], color = 'darkblue', linewidth=linewidth, label='Phage transcripts')
    ax.plot(cdf_list[1][0], cdf_list[1][1], color = 'lightblue', linewidth=linewidth, label='Non-phage transcripts')

    fig.legend(loc=(0.07, 0.86), frameon=False, fontsize=fontsize, ncol = 3)
    plt.savefig(os.path.join(args.outdir, "%s_at_cdf.%s"  %  (name, args.format)) , format = args.format)
    plt.clf()
    
for cdf_list, name in zip( [transcript_cdf, upstream_cdf], ['transcript', 'upstream'] ):
    draw_cdf(cdf_list, name, 30, 5)

        
    



