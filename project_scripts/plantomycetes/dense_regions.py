import argparse
import os
import sys
import numpy as np;
from collections import defaultdict
from pybedtools import BedTool, Interval
from Bio import SeqIO
import matplotlib.pyplot as plt;

from os import listdir

from afbio.sequencetools import get_at_content, sliding_window, array2fixed_length
from afbio.numerictools import find_elements_order





parser = argparse.ArgumentParser(description='Analyses AT profiles of the provided transcriptomes to predict those with silencers');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome folder");
parser.add_argument('--upstream', nargs = '?', default=50, type = int, help = "Upstream area to TSS");
parser.add_argument('--downstream', nargs = '?', default=20, type = int, help = "Downstream area to TSS");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
args = parser.parse_args();


P_RANGE = [1,3,5,7,10]
L_RANGE = len(P_RANGE)

def score_genome(mylist):
    return sum(mylist[1:1+L_RANGE]) + 0.4*sum(mylist[1+L_RANGE:1+2*L_RANGE])
    
 
def get_tss_at_contents(genes, fasta, upstream, downstream):
    #total_at = get_at_content("".join([str(x.seq) for x in fasta.values()]))
    tss = [fasta[x.chrom].seq[x.start-upstream: x.start+downstream] for x in genes if x.strand == '+']
    tss.extend([fasta[x.chrom].seq[x.end-downstream: x.end+upstream] for x in genes if x.strand == '-'])
    tss_at = np.array([get_at_content(x) for x in tss if x])
    
    median = np.median(tss_at);
    high_at = [np.percentile(tss_at, 100-p) - median for p in P_RANGE]
    low_at = [median - np.percentile(tss_at, p) for p in P_RANGE]
    
    return high_at, low_at;
    #print(len(tss_at))
    #print(high_at);
    #print(low_at);
    


gff_files = sorted([os.path.join(args.path, f) for f in listdir(args.path) if os.path.isfile(os.path.join(args.path, f)) and f.endswith('gff')])
fasta_files = sorted([os.path.join(args.path, f) for f in listdir(args.path) if os.path.isfile(os.path.join(args.path, f)) and f.endswith('fna')])


names = [];
high_at_list = [[] for i in P_RANGE];
#print(len(high_at_list))
differences_list = [[] for i in P_RANGE];


for gff, fna in zip(gff_files, fasta_files):
    names.append(os.path.basename(fna)[:-4])
    fasta = SeqIO.to_dict(SeqIO.parse(fna, "fasta"))
    genes = BedTool(gff)
    high_at, low_at = get_tss_at_contents(genes, fasta, args.upstream, args.downstream)

    




