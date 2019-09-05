import argparse
import os
import sys
import numpy as np;
from collections import defaultdict
from pybedtools import BedTool, Interval
from Bio import SeqIO
import matplotlib.pyplot as plt;

from afbio.sequencetools import get_at_content, sliding_window, array2fixed_length




parser = argparse.ArgumentParser(description='Explores relation between AT stretches and binding peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts coordinates, bed/gff file");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--length', nargs = '?', default=200, type = int, help = "Length of the upstream/downstream regions");
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
            

genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
transcripts = BedTool(args.transcripts);
phages = BedTool(args.phages);
regions = BedTool(args.path)
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
tr_intervals = [upstreams, transcripts, downstreams]

for ph_name, ptr in zip(ph_names, phaged_regions):
    for tr_name, intervals in zip(tr_names, tr_intervals):
        for region in ptr.intersect(b=intervals, f=0.5, u=True):
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


def get_transcript_profile(regions, transcripts, weighted=False):
    temp_dict = defaultdict(list);
    profiles = [];
    
    for region in regions.intersect(b=transcripts, f=0.5, wb=True):
        cds_name = dict([ x.split("=") for x in region[OFFSET+8].split(";") if x]) ['ID']
        temp_dict[cds_name].append(region);
        
    for temp_regions in temp_dict.values():
        local_profiles = []
        for region in temp_regions:
            if(weighted):
                count = float(region.attrs['topcoverage'])
            else:
                count = 1;
            tr_start = int(region[OFFSET+3])-1
            tr_end = int(region[OFFSET+4])
            tr_len = tr_end-tr_start          
            l_start = max(region.start - tr_start, 0)
            l_stop = min(region.end - tr_start, tr_len)
            arr = [0]*l_start + [count]*(l_stop-l_start) + [0]*(tr_len-l_stop)
            if(region[OFFSET+6] == '-'):
                arr = arr[::-1]
            local_profiles.append(array2fixed_length(arr, 100))
            
        local_profiles = np.array(local_profiles)
        profiles.append(local_profiles.mean(axis=0));
    
    profiles = np.array(profiles);
    return profiles.mean(axis=0);


def draw_transcript_profiles(transcript_profiles, mode, fontsize=28, linewidth=5):
    labels = ['phage', 'non-phage']
    colors = ['darkblue', 'lightblue']
    x_range = range(1, 101);

    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

    ax.set_xlabel("transcript 5\'->3\'[% of length]", fontsize=fontsize)
    if(mode=='plain'):
        ax.set_ylabel('Peak density', fontsize=fontsize)
    else:
        ax.set_ylabel('Peak coverage', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)
    
    for profile, color, label in zip(transcript_profiles, colors, labels):
        ax.plot(x_range, profile, color = color, linewidth=linewidth, label=label)

    fig.legend(loc=(0.15, 0.86), frameon=False, fontsize=fontsize, ncol = 2)
    plt.savefig(os.path.join(args.outdir, "transcript_profile_%s.%s"  %  (mode, args.format)) , format = args.format)
    plt.clf()


        
if(args.full):
    transcript_profiles_plain = [];
    transcript_profiles_weighted = [];
    for k, v in transcriptome_region_dict.items():
        if(k[1] == "transcript"):
            transcript_profiles_plain.append(get_transcript_profile(v, transcripts, weighted=False))
            transcript_profiles_weighted.append(get_transcript_profile(v, transcripts, weighted=True))

    for transcript_profiles, mode in zip([transcript_profiles_plain, transcript_profiles_weighted], ['plain', 'weighted']):
        draw_transcript_profiles(transcript_profiles, mode, fontsize=28, linewidth=5);
        
        
        
        
###################################################################################################
### upstream profile

def get_upstream_profile(regions, upstreams, weighted=False):
    profiles = [];
    
    for region in regions.intersect(b=upstreams, f=0.5, wb=True):
        if(weighted):
            count = float(region.attrs['topcoverage'])
        else:
            count = 1;
            
        tr_start = int(region[OFFSET+1])
        tr_end = int(region[OFFSET+2])
        tr_len = tr_end-tr_start          
        l_start = max(region.start - tr_start, 0)
        l_stop = min(region.end - tr_start, tr_len)
        
        arr = [0]*l_start + [count]*(l_stop-l_start) + [0]*(tr_len-l_stop)
        if(region[OFFSET+5] == '+'):
            arr = arr[::-1]
        profiles.append(arr);
        
    profiles = np.array(profiles);
    return profiles.mean(axis=0);



def draw_upstream_profiles(upstream_profiles, mode, fontsize=28, linewidth=5):
    labels = ['phage', 'non-phage']
    colors = ['darkblue', 'lightblue']
    x_range = range(len(upstream_profiles[0]));

    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

    ax.set_xlabel("distance from CDS [nt]", fontsize=fontsize)
    if(mode=='plain'):
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
    plt.savefig(os.path.join(args.outdir, "upstream_profile_%s.%s"  %  (mode, args.format)) , format = args.format) 
    plt.clf()






upstream_profiles_plain = [];
upstream_profiles_weighted = [];
for k, v in transcriptome_region_dict.items():
    if(k[1] == "upstream"):
        upstream_profiles_plain.append(get_upstream_profile(v, upstreams, weighted=False))
        upstream_profiles_weighted.append(get_upstream_profile(v, upstreams, weighted=True))

for upstream_profiles, mode in zip([upstream_profiles_plain, upstream_profiles_weighted], ['plain', 'weighted']):
    draw_upstream_profiles(upstream_profiles, mode, fontsize=28, linewidth=5);
        


        

        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


