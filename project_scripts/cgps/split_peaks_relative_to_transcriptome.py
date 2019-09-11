import argparse
import os
import sys
import numpy as np;
from collections import defaultdict
from pybedtools import BedTool, Interval
from Bio import SeqIO
import matplotlib.pyplot as plt;
import copy

from afbio.sequencetools import get_at_content, sliding_window, array2fixed_length, coverage2dict




parser = argparse.ArgumentParser(description='Explores positions of binding peaks on a genome');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks, gff format");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts coordinates, gff file");
parser.add_argument('--coverage', nargs = '?', required=True, type = str, help = "Path to the normalized coverage track, bed format");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--length', nargs = '?', default=200, type = int, help = "Length of the upstream/downstream regions");
parser.add_argument('--zscore', nargs = '?', default=2, type = int, help = "Z-score threshold for the binding peaks");
parser.add_argument('--at_length', nargs = '?', default=30, type = int, help = "Length of an AT-rich motif");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
parser.add_argument('--full', nargs = '?', default=False, const=True, type = bool, help = "If set, the full analyses is performed")
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
            
coverage = coverage2dict(args.coverage)
genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
transcripts = BedTool(args.transcripts);
phages = BedTool(args.phages);
regions = BedTool(args.path)
regions = BedTool([x for x in regions if float(x.score)>args.zscore])
OFFSET = regions.field_count()
up_downs = [get_upstream_downstream(x, genome, args.length) for x in transcripts]
upstreams = BedTool([x[0] for x in up_downs if x[0]])
downstreams = BedTool([x[1] for x in up_downs if x[1]])
phaged_regions = [regions.intersect(b=phages, u =True, f = 0.5), regions.intersect(b=phages, v =True, f = 0.5)]

print([len(x) for x in phaged_regions]);


already_discovered = set()
transcriptome_region_dict = defaultdict(list);
ph_names = ['phage', 'non-phage']
tr_names = ['upstream', 'transcript', 'downstream']
tr_fractions = [0.2, 0.2, 0.2]
tr_intervals = [upstreams, transcripts, downstreams]

for ph_name, ptr in zip(ph_names, phaged_regions):
    for tr_name, tr_fraction, intervals in zip(tr_names, tr_fractions, tr_intervals):
        for region in ptr.intersect(b=intervals, f=tr_fraction, u=True):
            if(region.name not in already_discovered):
                transcriptome_region_dict[(ph_name, tr_name)].append(region);
            already_discovered.add(region.name)
    for region in ptr:
        if(region.name not in already_discovered):
            transcriptome_region_dict[(ph_name, 'intergenic')].append(region);
    
  
for k, v in transcriptome_region_dict.items():
    transcriptome_region_dict[k] = BedTool(v)
    #print("\t".join(k), "\t", len(v))
        
    
        
###################################################################################################
### transcript profile


def get_transcript_profile(regions, transcripts, coverage, fraction):
    temp_dict = defaultdict(list);
    profiles = [];
    
    for transcript in transcripts.intersect(b=regions, F=fraction, u=True):
        temp_dict[transcript.name].append(transcript);
        
    for temp_transcripts in temp_dict.values():
        local_profiles = []
        for transcript in temp_transcripts:
            arr = coverage[transcript.chrom][transcript.start: transcript.stop]
            if(transcript.strand == '-'):
                arr = arr[::-1]
            local_profiles.append(array2fixed_length(arr, 100))
        local_profiles = np.array(local_profiles)
        profiles.append(local_profiles.mean(axis=0));
    
    profiles = np.array(profiles);
    return profiles.mean(axis=0);


def draw_transcript_profiles(transcript_profiles, fontsize=28, linewidth=5):
    labels = ['phage', 'non-phage']
    colors = ['darkblue', 'lightblue']
    x_range = range(1, 101);

    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

    ax.set_xlabel("transcript 5\'->3\'[% of length]", fontsize=fontsize)
    ax.set_ylabel('Peak coverage', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)
    
    for profile, color, label in zip(transcript_profiles, colors, labels):
        ax.plot(x_range, profile, color = color, linewidth=linewidth, label=label)

    fig.legend(loc=(0.15, 0.86), frameon=False, fontsize=fontsize, ncol = 2)
    plt.savefig(os.path.join(args.outdir, "transcript_profile.%s"  %  args.format) , format = args.format)
    plt.clf()


        
if(args.full):
    transcript_profiles = [];
    #transcript_profiles_weighted = [];
    for k, v in transcriptome_region_dict.items():
        if(k[1] == "transcript"):
            transcript_profiles.append(get_transcript_profile(v, transcripts, coverage, fraction=tr_fractions[1]))
            #transcript_profiles_weighted.append(get_transcript_profile(v, transcripts, weighted=True))
    draw_transcript_profiles(transcript_profiles, fontsize=28, linewidth=5);
        
        
        
        
###################################################################################################
### upstream profile

def get_upstream_profile(regions, upstreams, coverage, fraction, plain=False):
    profiles = [];
    for upstream in upstreams.intersect(b=regions, F=fraction, u=True):
        arr = coverage[upstream.chrom][upstream.start: upstream.end]
        if(plain):
            norma = max(arr)
            arr = [x/norma for x in arr]
        if(upstream.strand == '+'):
            arr = arr[::-1]
        profiles.append(copy.copy(arr));
        
    profiles = np.array(profiles);
    return profiles.mean(axis=0);



def draw_upstream_profiles(upstream_profiles, fontsize=28, linewidth=5, plain=False):
    labels = ['phage', 'non-phage']
    colors = ['darkblue', 'lightblue']
    x_range = range(len(upstream_profiles[0]));

    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

    ax.set_xlabel("distance from CDS [nt]", fontsize=fontsize)
    if(plain):
        ax.set_ylabel('Peak density', fontsize=fontsize)
    else:
        ax.set_ylabel('Peak coverage', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)
    
    for profile, color, label in zip(upstream_profiles, colors, labels):
        ax.plot(x_range, profile, color = color, linewidth=linewidth, label=label)

    fig.legend(loc=(0.15, 0.86), frameon=False, fontsize=fontsize, ncol = 2)
    if(plain):
        plt.savefig(os.path.join(args.outdir, "upstream_profile_plain.%s"  %  args.format) , format = args.format)
    else:
        plt.savefig(os.path.join(args.outdir, "upstream_profile.%s"  %  args.format) , format = args.format)
    plt.clf()





if(args.full):
    upstream_profiles = [];
    upstream_profiles_plain = [];
    for k, v in transcriptome_region_dict.items():
        if(k[1] == "upstream"):
            upstream_profiles.append(get_upstream_profile(v, upstreams, coverage, fraction=tr_fractions[0], plain=False))
            upstream_profiles_plain.append(get_upstream_profile(v, upstreams, coverage, fraction=tr_fractions[0], plain=True))


    draw_upstream_profiles(upstream_profiles, fontsize=28, linewidth=5, plain=False);
    draw_upstream_profiles(upstream_profiles_plain, fontsize=28, linewidth=5, plain=True)
    
        
        
        
###################################################################################################
### AT at_content and topcoverage

def get_peak_at_content_exact(regions, genome, at_length):
    tail = at_length//2
    res = [];
    for region in regions:
        chrom = genome[region.chrom];
        middle = int(region.name)
        seq = str(chrom[middle-tail:middle+tail].seq.upper())
        res.append(get_at_content(seq));
    return res;

def get_peak_at_content_max(regions, genome, at_length):
    res = [];
    for region in regions:
        chrom = genome[region.chrom];
        seq = str(chrom[region.start:region.stop].seq.upper())
        max_at = max([get_at_content(x) for x in sliding_window(seq, at_length)])
        res.append(max_at);
    return res;


def get_peak_topcoverage(regions):
    return [float(x.attrs['topcoverage']) for x in regions];


def draw_barplot(data, name, ylabel, labels, fontsize=28, linewidth = 5, width = 0.25, ylim=False):
    
    arrsize = len(data[0])
    means = [[],[]]
    yerr_list = [[],[]]
    
    for i, list2d in enumerate(data):
        for list1d in list2d:
            #print(len(list1d))
            mean = np.mean(list1d)
            means[i].append(mean)
            yerr_list[i].append(( mean - np.percentile(list1d, 25), np.percentile(list1d, 75) - mean ))
    yerr_list = [np.array(x).transpose() for x in yerr_list]

    x = np.arange(arrsize)
    barlabels = ['phage', 'non-phage']
    barcolors = ['darkblue', 'lightblue']

    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.95])

    ax.set_xlabel('Type of target', fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)
    if(ylim):
        bottom = min([min(x) for x in means])
        ax.set_ylim(bottom=bottom);
    
    for adj, yvals, yerr, label, color in zip( [-width/2, width/2], means, yerr_list, barlabels, barcolors):
        #print(yerr)
        ax.bar(x + adj, yvals, width, label=label, color = color)
        plt.errorbar(x + adj, yvals, yerr=yerr, linestyle='', capsize=12, linewidth=linewidth/2, capthick=linewidth/2, color='black')
    #print()
    
    fig.legend(loc=(0.2, 0.9), frameon=False, fontsize=fontsize, ncol = 2)
    if(ylim):
        plt.savefig(os.path.join(args.outdir, "peaks_%s_ylim.%s"  % (name, args.format)) , format = args.format)
    else:
        plt.savefig(os.path.join(args.outdir, "peaks_%s.%s"  % (name, args.format)) , format = args.format)
    plt.clf()





at_content_exact = []
topcoverage_list = []
at_content_max = []
tr_names = ['upstream', 'transcript', 'downstream', 'intergenic']

for ph_name in ph_names:
    at_content_exact.append([get_peak_at_content_exact(transcriptome_region_dict[(ph_name, x)], genome, args.at_length) for x in tr_names])
    at_content_max.append([get_peak_at_content_max(transcriptome_region_dict[(ph_name, x)], genome, args.at_length) for x in tr_names])
    topcoverage_list.append([get_peak_topcoverage(transcriptome_region_dict[(ph_name, x)]) for x in tr_names])
    
draw_barplot(at_content_exact, "at_exact", "AT content", tr_names, fontsize=28, linewidth = 5, width = 0.25, ylim=False)
draw_barplot(at_content_max, "at_max", "AT content", tr_names, fontsize=28, linewidth = 5, width = 0.25, ylim=False)
draw_barplot(topcoverage_list, "topcoverage", "Top coverage [nomalized]", tr_names, fontsize=28, linewidth = 5, width = 0.25, ylim=False)

        


        

        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


