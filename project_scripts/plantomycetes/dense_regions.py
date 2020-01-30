import argparse
import os
import sys
import numpy as np;
from collections import defaultdict
from pybedtools import BedTool, Interval
from Bio import SeqIO
import matplotlib.pyplot as plt;

from os import listdir

from afbio.sequencetools import get_at_content
from afbio.numerictools import find_elements_order





parser = argparse.ArgumentParser(description='Analyses AT profiles of the provided transcriptomes to predict those with silencers');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome folder");
parser.add_argument('--score_table', nargs = '?', required=True, type = str, help = "Path to the score table")
parser.add_argument('--upstream', nargs = '?', default=50, type = int, help = "Upstream area to TSS");
parser.add_argument('--downstream', nargs = '?', default=20, type = int, help = "Downstream area to TSS");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
args = parser.parse_args();


def get_at_smooth_profile(fasta, length=100):
    profile = []
    step = length//2
    for seqrecord in fasta.values():
        for start in range(0, len(seqrecord)-step, step):
            stop = start + length
            seq = str(seqrecord.seq)
            profile.append(get_at_content(seq))
    norma = np.mean(profile)
    return np.array(profile)/norma
    

def get_density_score(profile):
    temp = [float(x) for x in profile]
    return sum(temp)/len(temp)
    
 
def get_profile(genes, fasta, upstream, downstream):
    #total_at = get_at_content("".join([str(x.seq) for x in fasta.values()]))
    tss = [ (x.chrom, x.start-upstream, fasta[x.chrom].seq[x.start-upstream: x.start+downstream]) for x in genes if x.strand == '+']
    tss.extend( [(x.chrom, x.end-downstream, fasta[x.chrom].seq[x.end-downstream: x.end+upstream]) for x in genes if x.strand == '-'])
    tss_start_at = [ (x[0], x[1], get_at_content(x[2])) for x in tss if x[2] ]
    tss_at  = [x[2] for x in tss_start_at]
    threshold = np.percentile(tss_at, 90)
    total_tss = [ (x[0], x[1]) for x in tss_start_at]
    passed_tss = [ (x[0], x[1]) for x in tss_start_at if x[2]>=threshold]
    
    length = 100000
    step = length//2
    profile = []
    for chrom, seqrecord in fasta.items():
        #print(len(seqrecord));
        for start in range(0, len(seqrecord)-step, step):
            stop = start + length
            local_passed = 0;
            for local_chrom, local_start in passed_tss:
                if(chrom == local_chrom and local_start>=start and local_start<stop):
                    local_passed += 1;
            total_passed = 1;
            for local_chrom, local_start in total_tss:
                if(chrom == local_chrom and local_start>=start and local_start<stop):
                    total_passed += 1;
                    
            profile.append(local_passed/total_passed)
    norma = np.mean(profile)
    profile = np.array(profile)/norma
    
    return tuple(list(sorted(profile, reverse = True))[:10]), len(profile), profile
    
            #print(start, stop, local_passed);
    
    
    
    #print(len(passed_tss), len(genes))

    


gff_files = sorted([os.path.join(args.path, f) for f in listdir(args.path) if os.path.isfile(os.path.join(args.path, f)) and f.endswith('gff')])
fasta_files = sorted([os.path.join(args.path, f) for f in listdir(args.path) if os.path.isfile(os.path.join(args.path, f)) and f.endswith('fna')])



name2profile = {};
name2smooth = {};
for gff, fna in zip(gff_files, fasta_files):
    name = os.path.basename(fna)[:-4]
    fasta = SeqIO.to_dict(SeqIO.parse(fna, "fasta"))
    genes = BedTool(gff)
    profile, lp, raw = get_profile(genes, fasta, args.upstream, args.downstream)
    name2profile[name] = ["%1.2f" % x for x in profile], lp, raw
    #name2smooth[name] = get_at_smooth_profile(fasta, length=400)
    #print(["%1.2f" % x for x in profile])
    sys.stderr.write("%s\n" % name)
 

final_table = []
with open(args.score_table) as f:
    next(f)
    for l in f:
        a = l.strip().split("\t")
        name = a[0]
        tss_score = float(a[1])
        profile, lp, raw = name2profile[name]
        density_score = get_density_score(profile)
        score = tss_score*(1+0.75*density_score)
        final_table.append((name, score, tss_score, density_score))
        
final_table.sort(key = lambda x: x[1], reverse = True);

print("\t".join(("Name", "Total Score", "AT Score", "Density Score")))
for mylist in final_table:
    print( "%s\t%1.1f\t%1.2f\t%1.2f" % mylist)
    
    
    

    
############################# DRAWING SECTION #############################
def draw(bars, name, width = 1, fontsize = 24, linewidth = 4):



    fig, ax = plt.subplots(figsize = (16, 9))
    plt.tight_layout(rect=[0.08, 0.08, 1, 1])
    x = np.arange(len(bars))
    ax.bar(x, bars, width, color = 'lightblue')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('0->100% genome length', fontsize=fontsize)
    ax.set_ylabel('number of AT-rich genes (normalized)', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize-4, top=False, right=False, bottom=False)
    plt.xticks([])

    #plt.show()
    plt.savefig(os.path.join(args.outdir, "at_profile.%s.%s") % (name, args.format) , format = args.format)
    plt.close()
    
    
#for mylist in final_table[:5]:
    #name = mylist[0]
    #profile = name2smooth[name]
    #draw(profile, name)


