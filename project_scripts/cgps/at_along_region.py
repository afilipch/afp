import argparse
import os
import sys
import numpy as np;
from collections import defaultdict
from pybedtools import BedTool, Interval
from Bio import SeqIO
import matplotlib.pyplot as plt;

from afbio.sequencetools import get_at_content, sliding_window, array2fixed_length, coverage2dict

#from afbio.sequencetools import get_at_content, sliding_window, array2fixed_length, coverage2dict




parser = argparse.ArgumentParser(description='Checks the AT content along the regions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks, gff format");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
parser.add_argument('--zscore', nargs = '?', default=2, type = int, help = "Z-score threshold for the binding peaks");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
args = parser.parse_args();


def get_at_profile(regions, genome):
    res = [];
    for region in regions:
        chrom = genome[region.chrom];
        seq = str(chrom[region.start-10:region.stop+10].seq.upper())
        res.append([get_at_content(x) for x in sliding_window(seq, 20)])
    
    res = np.array(res)
    return res.mean(axis=0)



genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
regions = BedTool(args.path)
regions = BedTool([x for x in regions if float(x.score)>args.zscore])
phages = BedTool(args.phages);

phaged_regions = [regions.intersect(b=phages, u =True, f = 0.5), regions.intersect(b=phages, v =True, f = 0.5)]
profiles = []
for ptr in phaged_regions:
    profiles.append(get_at_profile(ptr, genome))
    
    
def draw_profiles(profiles, fontsize=28, linewidth=5):
    labels = ['phage', 'non-phage']
    colors = ['darkblue', 'lightblue']
    x_range = range(len(profiles[0]));

    fig, ax = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

    ax.set_xlabel("position centered to peak [nt]", fontsize=fontsize)
    ax.set_ylabel('AT content', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)
    
    for profile, color, label in zip(profiles, colors, labels):
        ax.plot(x_range, profile, color = color, linewidth=linewidth, label=label)

    fig.legend(loc=(0.15, 0.86), frameon=False, fontsize=fontsize, ncol = 2)
    plt.savefig(os.path.join(args.outdir, "regions_at_profile.%s"  %  args.format) , format = args.format) 
    plt.clf()
    
    
draw_profiles(profiles, fontsize=28, linewidth=5)    
    

