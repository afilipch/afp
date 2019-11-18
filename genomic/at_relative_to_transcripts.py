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

from afbio.sequencetools import get_at_content, sliding_window, array2fixed_length
from afbio.numerictools import CDF




parser = argparse.ArgumentParser(description='Analyses AT content relative to the transcripts and phages positions on a genome');
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts coordinates, bed/gff file");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--window', nargs = '?', default=40, type = int, help = "Local length of the widonwt to check for an AT content");
parser.add_argument('--lookup', nargs = '?', default=500, type = int, help = "Total length of the segment upstream a transcript to check for an AT content");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");

#parser.add_argument('--tss', nargs = '?', required=True, type = str, help = "Path to the TSS annotation file");
#parser.add_argument('--stranded', nargs = '?', const=True, default=False, type = bool, help = "If set, the data considered to be stranded");
args = parser.parse_args();


def transcript_content(interval, genome, window):
    #print(type(interval.start))
    if(interval.strand == '+'):
        seq = str(genome[interval.chrom][interval.start:interval.stop].seq.upper())
    elif(interval.strand == '-'):
        seq = str(genome[interval.chrom][interval.start:interval.stop].seq.reverse_complement().upper())
        
    profile = [get_at_content(x) for x in sliding_window(seq, window)]
    return array2fixed_length(profile, 100)
        
    #return profile


def transcript2upstream(interval, genome, window, lookup):
    
    seq = None
    chrom = genome[interval.chrom]
    if(interval.strand == '+'):
        start = interval.start - lookup - window;
        end = interval.start;
        if(start >= 0 ):
            seq = str(chrom[start:end].seq.reverse_complement().upper())
    elif(interval.strand == '-'):
        start = interval.stop
        end = interval.stop + lookup + window;
        if(end<len(chrom)):
            seq = str(chrom[start:end].seq.upper())
    if(seq):
        return [get_at_content(x) for x in sliding_window(seq, window)]
    




genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
transcripts = BedTool([x for x in BedTool(args.transcripts) if len(x) > 100]);
phages = BedTool(args.phages);


phaged_transcripts = [transcripts.intersect(b=phages, u =True, f = 0.5), transcripts.intersect(b=phages, v =True, f = 0.5)]
transcript_profiles = []
upstream_profiles = []
transcript_at_content = []
upstream_at_content = []

for ptr in phaged_transcripts:
    upstream_local_profiles = [transcript2upstream(x, genome, args.window, args.lookup) for x in ptr]
    upstream_local_profiles = [x for x in upstream_local_profiles if x]
    upstream_at_content.append([x[0] for x in upstream_local_profiles]);
    upstream_local_profiles = np.array(upstream_local_profiles)
    upstream_profiles.append(upstream_local_profiles.mean(axis=0))
    

    transcript_local_profiles = []
    temp_dict = defaultdict(list);
    for transcript in ptr:
        #print(transcript)
        temp_dict[transcript.name].append(transcript_content(transcript, genome, args.window))
        
    for l in temp_dict.values():
        arr = np.array(l);
        transcript_local_profiles.append(arr.mean(axis=0))
    
    transcript_local_profiles = np.array(transcript_local_profiles)
    transcript_profiles.append(transcript_local_profiles.mean(axis=0));
    transcript_at_content.append(transcript_local_profiles.mean(axis=1));
    
print([len(x) for x in transcript_at_content])
upstream_cdf = [CDF(x) for x in upstream_at_content]
transcript_cdf = [CDF(x) for x in transcript_at_content]
        
        
        
    #at_content = [np.mean(x) for x in temp_dict.values()]
    #transcript_cdf.append(CDF(at_content, zerovalue=0))
    #at_content_list.append(np.array(at_content))
    
    


#upstream_cdf = [];
#for ptr in phaged_transcripts:
    #at_content = [transcript_upstream(x, genome, 20, 60) for x in ptr];
    #at_content = [max(x) for x in at_content if x]
    #upstream_cdf.append(CDF(at_content, zerovalue=0))
    #at_content_list.append(np.array(at_content))
    
#print("type of regions", "num of regions", "mean AT[%]", "median AT[%]", "std AT[%]")
#for at_content, name in zip(at_content_list, ["phage transcripts", "non-phage transcripts", "phage upstream", "non-phage upstream"]):
    #print("%s\t%d\t%1.1f\t%1.1f\t%1.1f" % (name, len(at_content), at_content.mean()*100, np.median(at_content)*100, at_content.std()*100 ))
    
#print(upstream_cdf[1])
#print(upstream_cdf[0])


############################################################################################################
### Plotting section

### Plot phage versus transcript distribution

def draw_cdf(cdf_list, name, fontsize=28, linewidth=5):
    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

    ax.set_xlabel('AT-content', fontsize=fontsize)
    ax.set_ylabel('CDF', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)

    ax.plot(cdf_list[0][0], cdf_list[0][1], color = 'darkblue', linewidth=linewidth, label='phage')
    ax.plot(cdf_list[1][0], cdf_list[1][1], color = 'lightblue', linewidth=linewidth, label='non-phage')

    fig.legend(loc=(0.15, 0.86), frameon=False, fontsize=fontsize, ncol = 2)
    plt.savefig(os.path.join(args.outdir, "%s_at_cdf.%s"  %  (name, args.format)) , format = args.format)
    plt.clf()
    
    
def draw_transcript_profiles(transcript_profiles, fontsize=28, linewidth=5):
    labels = ['phage', 'non-phage']
    colors = ['darkblue', 'lightblue']
    x_range = range(1, 101);

    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

    ax.set_xlabel("transcript 5\'->3\'[% of length]", fontsize=fontsize)
    ax.set_ylabel('AT-content', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)
    
    for profile, color, label in zip(transcript_profiles, colors, labels):
        ax.plot(x_range, profile, color = color, linewidth=linewidth, label=label)

    fig.legend(loc=(0.15, 0.86), frameon=False, fontsize=fontsize, ncol = 2)
    plt.savefig(os.path.join(args.outdir, "transcript_at_profile.%s"  %  args.format) , format = args.format)
    plt.clf()
    
    
def draw_upstream_profiles(upstream_profiles, fontsize=28, linewidth=5):
    labels = ['phage', 'non-phage']
    colors = ['darkblue', 'lightblue']
    x_range = range(len(upstream_profiles[0]));

    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

    ax.set_xlabel("distance from CDS [nt]", fontsize=fontsize)
    ax.set_ylabel('AT-content', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)
    
    for profile, color, label in zip(upstream_profiles, colors, labels):
        ax.plot(x_range, profile, color = color, linewidth=linewidth, label=label)

    fig.legend(loc=(0.15, 0.86), frameon=False, fontsize=fontsize, ncol = 2)
    plt.savefig(os.path.join(args.outdir, "upstream_at_profile.%s"  %  args.format) , format = args.format)
    plt.clf()
    
#for cdf_list, name in zip( [transcript_cdf, upstream_cdf], ['transcript', 'upstream'] ):
    #draw_cdf(cdf_list, name, 30, 5)
    
    
draw_transcript_profiles(transcript_profiles, fontsize=28, linewidth=5)
draw_upstream_profiles(upstream_profiles, fontsize=28, linewidth=5)
draw_cdf(upstream_cdf, 'upstream', fontsize=28, linewidth=5)
draw_cdf(transcript_cdf, 'transcript', fontsize=28, linewidth=5)

        
    



