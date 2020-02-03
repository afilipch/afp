#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores relation between GC drops and binding peaks'''

import argparse
import os
import sys
import numpy as np;
from scipy.stats import spearmanr, pearsonr
from collections import defaultdict;
import pandas as pd;
from pybedtools import BedTool, Interval
from Bio import SeqIO

from afbio.sequencetools import sliding_window


parser = argparse.ArgumentParser(description='Explores relation between GC drops and binding peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks");
parser.add_argument('--gcdrops', nargs = '?', required=True, type = str, help = "Path to the GC drops");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format");
parser.add_argument('--length', nargs = '?', default=7, type = int, help = "Length of the explored binding motifs");
parser.add_argument('--minz', nargs = '?', default=3, type = float, help = "Minimum required Z-score for the binding peaks");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory to store sequences of bound and unbound GC drops");

args = parser.parse_args();


genome = next(SeqIO.parse(args.genome, 'fasta')).seq.upper();
    
def count_nmers(sequences, length):
    counter = defaultdict(int);
    for seq in sequences:
        for nmer in set(sliding_window(seq, length)):
            counter[nmer] += 1;
    return sorted(counter.items(), key = lambda x: x[1], reverse = True);
    

regions = BedTool(args.path)
sys.stderr.write("total binding regions: %d\n" % len(regions));
regions = BedTool([ x for x in regions if float(x.score) >= args.minz])
sys.stderr.write("filtered binding regions: %d\n\n" % len(regions));
gcdrops_intervals = BedTool(args.gcdrops);


gcdrops_regions = gcdrops_intervals.intersect(regions, u = True)
gcdrops_unbound = gcdrops_intervals.intersect(regions, v = True)
sys.stderr.write("total drops: %d\nbound drops: %d\n\n" % (len(gcdrops_intervals) ,len(gcdrops_regions)));

gcdrop2seq = dict( [(">%d_%d" % (x.start, x.end), str(genome[x.start: x.end])) for x in gcdrops_regions] )
with open(os.path.join(args.outdir, "bound_minZ_%1.1f.fasta" % args.minz), 'w') as f:
    for nameseq in gcdrop2seq.items():
        if(len(nameseq[1]) >= 8): 
            f.write("%s\n%s\n" % nameseq);

unbound_gcdrop2seq = dict( [(">%d_%d" % (x.start, x.end), str(genome[x.start: x.end])) for x in gcdrops_unbound] )
with open(os.path.join(args.outdir, "unbound_minZ_%1.1f.fasta" % args.minz), 'w') as f:
    for nameseq in unbound_gcdrop2seq.items():
        if(len(nameseq[1]) >= 8): 
            f.write("%s\n%s\n" % nameseq);       



nmercounts = count_nmers(gcdrop2seq.values(), args.length)
#sys.stderr.write("Number of gcdrops inside the binding regions\t%d\n\n" % len(gcdrops_regions))
for nmercount in nmercounts[:10]:
    sys.stderr.write("%s\t%d\n" % nmercount)








































