import argparse
import os
import sys
import numpy as np;
from collections import defaultdict
from pybedtools import BedTool, Interval
from Bio import SeqIO
import matplotlib.pyplot as plt;

import numpy as np

from afbio.generators import get_only_files
from afbio.at_stretches import get_at_rich_stretches, stretch_score



parser = argparse.ArgumentParser(description='Extracts top AT-rich regions around TSS');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome folder");
parser.add_argument('--upstream', nargs = '?', default=75, type = int, help = "Upstream area to TSS");
parser.add_argument('--downstream', nargs = '?', default=25, type = int, help = "Downstream area to TSS");
parser.add_argument('--top', nargs = '?', default=10, type = int, help = "Output top N stretches");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
args = parser.parse_args();


 
def get_tss_at_contents(genes, fasta, upstream, downstream, top):
    stretches = []
    tss_list = [ (x, fasta[x.chrom].seq[x.start-upstream: x.start+downstream]) for x in genes if x.strand == '+']
    tss_list.extend([ (x, fasta[x.chrom].seq[x.end-downstream: x.end+upstream]) for x in genes if x.strand == '-' ])
    for tss, seq in tss_list:
        for adstart, adend, lseq, at_fraction, gc_count in get_at_rich_stretches(seq, 5, 20, 4, 0.4):
            stretches.append(( tss.chrom, tss.strand, tss.start + adstart, tss.start + adend, lseq, stretch_score(at_fraction, adend, adstart), seq )) 
    
    stretches.sort(key = lambda x: x[5], reverse = True)
    return stretches[:top]

    
raw_files = sorted(get_only_files(args.path))
files_list = [ (os.path.basename(x[1][0]).split(".")[0], x[1][0], x[1][1]) for x in enumerate(zip(raw_files, raw_files[1:])) if x[0] % 2 == 0]
#sys.exit()

for name, fasta_path, gff_path in files_list[:]:
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    genes = BedTool(gff_path)
    with open(os.path.join(args.outdir, "%s.at_stretches.tsv" % name), 'w') as f:
        f.write( "chromosome\tstrand\tstart\tstop\tseq\tscore\tall_tss_seq\n")
        for el in get_tss_at_contents(genes, genome, args.upstream, args.downstream, args.top):
            f.write( "%s\t%s\t%d\t%d\t%s\t%1.2f\t%s\n" % el)

    
        #for c, stretch in enumerate(stretches, start = 1):
            #f2.write( "%s\t%d\t%d\tstretch_%d\t%1.1f\t+\n" % (stretch[0], stretch[1], stretch[2], c, stretch[-1]))



