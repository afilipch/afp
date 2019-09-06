import argparse
import os
import sys
import numpy as np;
from collections import defaultdict
from pybedtools import BedTool, Interval
from Bio import SeqIO
import matplotlib.pyplot as plt;

#from afbio.sequencetools import get_at_content, sliding_window, array2fixed_length, coverage2dict




parser = argparse.ArgumentParser(description='Compares real peaks with random ones');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks, bed/gff format");
parser.add_argument('--control', nargs = '?', required=True, type = str, help = "Path to the control/random peaks, bed file");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts coordinates, gff file");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--length', nargs = '?', default=200, type = int, help = "Length of the upstream/downstream regions");
parser.add_argument('--zscore', nargs = '?', default=2, type = int, help = "Z-score threshold for the binding peaks");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
args = parser.parse_args();



#########################################################################################################################
### Data split section

def get_upstream_downstream(interval, genome, length):
    i1, i2 = None, None
    chrom = genome[interval.chrom]
    
    start = interval.start - length;
    stop = interval.start;
    if(start >= 0 ):
        i1 = Interval(interval.chrom, start, stop, interval.name, interval.score, interval.strand)
        
    start = interval.stop
    stop = interval.stop + length
    if(stop<len(chrom)):
        i2 = Interval(interval.chrom, start, stop, interval.name, interval.score, interval.strand)
        
    if(interval.strand == '-'):
        return i2, i1
    else:
        return i1, i2
    
    
    
genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
transcripts = BedTool(args.transcripts);
phages = BedTool(args.phages);
regions = BedTool(args.path)
regions = BedTool([x for x in regions if float(x.score)>args.zscore])
control = BedTool(args.control);

up_downs = [get_upstream_downstream(x, genome, args.length) for x in transcripts]
upstreams = BedTool([x[0] for x in up_downs if x[0]])
downstreams = BedTool([x[1] for x in up_downs if x[1]])

phaged_regions_list = [[regions.intersect(b=phages, u =True, f = 0.5), regions.intersect(b=phages, v =True, f = 0.5)], [control.intersect(b=phages, u =True, f = 0.5), control.intersect(b=phages, v =True, f = 0.5)]]


already_discovered = set()
transcriptome_region_dict = defaultdict(list);
r_names = ['real', 'control']
ph_names = ['phage', 'non-phage']
tr_names = ['upstream', 'transcript', 'downstream']
tr_intervals = [upstreams, transcripts, downstreams]

for r_name, phaged_regions in zip(r_names, phaged_regions_list):
    for ph_name, ptr in zip(ph_names, phaged_regions):
        for tr_name, intervals in zip(tr_names, tr_intervals):
            for region in ptr.intersect(b=intervals, f=0.5, u=True):
                if(region.name not in already_discovered):
                    transcriptome_region_dict[(r_name, ph_name, tr_name)].append(region);
                already_discovered.add(region.name)
        for region in ptr:
            if(region.name not in already_discovered):
                transcriptome_region_dict[(r_name, ph_name, 'intergenic')].append(region);

tr_names = ['upstream', 'transcript', 'downstream', 'intergenic']
for k, v in transcriptome_region_dict.items():
    print(k, len(v));
    
    
###get fraction of regions per transcript type;

fractions = [[]]*4
normas = len(regions), len(control)
for i, (norma, r_name) in enumerate(zip(normas, r_names)):
    for j, ph_name in enumerate(ph_names):
        fractions[i*2+j] = [len(transcriptome_region_dict[(r_name, ph_name, x)])/norma for x in tr_names]
    

###Plotting all data in one plot
    
def draw_barplot_all(data, labels, name, fontsize=28, linewidth = 5, width = 0.25):
    if(len(data) == 4):
         adjustmens = [-width*1.5, -width/2, width/2, width*1.5]
    elif(len(data) == 2):
        adjustmens = [-width/2, width/2]
            
        
    

    x = np.arange(len(data[0]))
    barcolors = ['darkblue', 'lightblue', 'darkgray', 'lightgray']
    barlabels = ['phage', 'non-phage', 'control phage', 'control non-phage']


    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.9])

    ax.set_xlabel('type of target', fontsize=fontsize)
    ax.set_ylabel('fraction of regions', fontsize=fontsize)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)    
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)

    
    for adj, yvals, label, color in zip( adjustmens, data, barlabels, barcolors):
        ax.bar(x + adj, yvals, width, label=label, color = color)

    
    fig.legend(loc=(0.25, 0.85), frameon=False, fontsize=fontsize, ncol = 2)
    plt.savefig(os.path.join(args.outdir, "peaks_vs_control_targets_%s.%s"  % (name, args.format)) , format = args.format)
    plt.clf()
    
    
draw_barplot_all(fractions, tr_names, 'all', fontsize=28, linewidth = 5, width = 0.16)



###Plotting phage/non-phage data in one plot

phage_fractions = [sum(x) for x in fractions[:2]], [sum(x) for x in fractions[2:]]

def draw_barplot_phage(data, fontsize=28, linewidth = 5, width = 0.25):

    x = np.arange(len(data[0]))
    barcolors = ['blue', 'gray']
    barlabels = ['regions', 'control']
    labels = ['phage', 'non-phage']


    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.95])

    ax.set_xlabel('type of target', fontsize=fontsize)
    ax.set_ylabel('fraction of regions', fontsize=fontsize)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)    
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(top=1);
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)

    
    for adj, yvals, label, color in zip( [-width/2, width/2], data, barlabels, barcolors):
        ax.bar(x + adj, yvals, width, label=label, color = color)

    
    fig.legend(loc=(0.25, 0.85), frameon=False, fontsize=fontsize, ncol = 2)
    plt.savefig(os.path.join(args.outdir, "peaks_vs_control_targets_phages.%s"  % args.format) , format = args.format)
    plt.clf()
    
draw_barplot_phage(phage_fractions, fontsize=28, linewidth = 5, width = 0.25)


###Plotting targets only for real data in one plot

regions_fractions = [x/sum(fractions[0]) for x in fractions[0]], [x/sum(fractions[1]) for x in fractions[1]]

print(regions_fractions)

draw_barplot_all(regions_fractions, tr_names, 'targets', fontsize=28, linewidth = 5, width = 0.16)
 
    






    






























