#! /usr/bin/python
'''Analyses TSS differences on phage versus non-phage'''

import argparse
import sys
import os
from collections import defaultdict, Counter
import matplotlib.pyplot as plt;
import numpy as np;
import random;

from Bio import SeqIO;
from pybedtools import BedTool, Interval

from afbio.sequencetools import get_at_content
from afbio.numerictools import CDF




parser = argparse.ArgumentParser(description='Analyses TSS differences on phage versus non-phage');
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts coordinates augmented with tss information [output of ~/afp/genomic/at_relative_to_transcripts.py], bed/gff file");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--length', nargs = '?',   default= 30, type = int, help = "Length of the segment right upstream");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();

transcripts = BedTool(args.transcripts);
phages = BedTool(args.phages);
phaged_transcripts = [transcripts.intersect(b=phages, u =True, f = 0.5), transcripts.intersect(b=phages, v =True, f = 0.5)]


############################################################################################################
### TSS count section


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


############################################################################################################
### AT-section

def get_control_transcripts(genome, distances):
    control = []
    chroms = list(genome.keys())
    tss_distance = max(distances)
    for variants in range(1, 5):
        for num in range(500):
            chrom = random.choice(chroms)
            length = len(genome[chrom])
            global_start = random.randint(tss_distance, length)
            end = global_start + 100
            for _ in range(variants):
                start = global_start - random.choice(distances)
                #print(chrom, start, start+100, "%d_%d" % (variants, num), '0', '+')
                control.append(Interval(chrom, start, end, "%d_%d" % (variants, num), '0', '+'));
    return control
                
                
                
            
    

def transcript2upstream(interval, genome, length):
    
    seq = None
    if(interval.strand == '+'):
        start = interval.start - length;
        end = interval.start;
        if(start >= 0 ):
            seq = str(genome[interval.chrom][start:end].seq.reverse_complement().upper())
    elif(interval.strand == '-'):
        start = interval.stop
        end = interval.stop + length;
        chrom = genome[interval.chrom]
        if(end<len(chrom)):
            seq = str(chrom[start:end].seq.upper())
    if(seq):
        return get_at_content(seq);
    else:
        return None
    

genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))


distances = [];
for ptr in phaged_transcripts:
    name2diff = defaultdict(list)
    for transcript in ptr:
        name2diff[transcript.name].append(transcript.start)
    distances.extend([ max(x) - min(x) for x in name2diff.values() if len(x) == 2]);
    

control = get_control_transcripts(genome, distances);
phaged_transcripts.append(control)


    
tss_at_variations = [];
tss_at_means = [];
tss_at_maxes = [];
tss_at_mins = [];
for ptr in phaged_transcripts:
    temp_dict = defaultdict(list);
    variants2variations = defaultdict(list)
    variants2means = defaultdict(list)
    variants2maxes = defaultdict(list)
    variants2mins = defaultdict(list)
    
    for transcript in ptr:
        at_content = transcript2upstream(transcript, genome, args.length)
        if(at_content):
            temp_dict[transcript.name].append(at_content)
            
            
    for at_contents in temp_dict.values():
        variants2variations[len(at_contents)].append(max(at_contents) - min(at_contents))
        variants2means[len(at_contents)].append(np.mean(at_contents))
        variants2maxes[len(at_contents)].append(max(at_contents))
        variants2mins[len(at_contents)].append(min(at_contents))
    tss_at_variations.append(variants2variations);
    tss_at_means.append(variants2means);
    tss_at_maxes.append(variants2maxes);
    tss_at_mins.append(variants2mins);

    
print("\nVariation\ntranscripts per gene\tnum genes phage\tnum genes non-phage\tnum genes random\tmean AT phage\tmean AT non-phage\tmean AT random\tvariation phage\tvariation non-phage\tvariation random\tmax AT phage\tmax AT non-phage\tmax AT random\tmin AT phage\tmin AT non-phage\tmin AT random");
for i in range(1, 5):
    a = tuple([i] + [len(x[i]) for x in tss_at_variations] + [np.mean(x[i]) for x in tss_at_means] + [np.mean(x[i]) for x in tss_at_variations] + [np.mean(x[i]) for x in tss_at_maxes] + [np.mean(x[i]) for x in tss_at_mins])
    print("%d\t%d\t%d\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f" %  a)
    

#print("\nVariation\ntranscripts per gene\tnum genes phage\tnum genes non-phage\tnum genes random\tmean AT phage\tmean AT non-phage\tmean AT random");    
#for i in range(1, 5):
    #a = tuple([i] + [len(x[i]) for x in tss_at_means] + [np.mean(x[i]) for x in tss_at_means])
    #print("%d\t%d\t%d\t%d\t%1.2f\t%1.2f\t%1.2f" %  a)
        


















    
