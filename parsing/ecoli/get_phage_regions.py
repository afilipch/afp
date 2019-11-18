'''Creates phage regions out of phage genes list and ncbi annotation'''

import argparse
import sys
import os
from pybedtools import BedTool



parser = argparse.ArgumentParser(description='Creates phage regions out of phage genes list and ncbi annotation');
parser.add_argument('--phage', nargs = '?', required=True, type = str, help = "Path to the phage genes file, custom format");
parser.add_argument('--ncbi', nargs = '?', required=True, type = str, help = "Path to the ncbi annotation");
parser.add_argument('--flank', nargs = '?', default=600, type = int, help = "flanks size for the phages on their 5prime end")
args = parser.parse_args();

genes = set();
tss_list = []
with open(args.phage) as f:
    for l in f:
        a = l.strip().split(" ")
        genes.add(a[0]);
        
cds_list = [x for x in BedTool(args.ncbi) if x[2] == 'CDS']
found = 0
for cds in cds_list:
    if(cds.attrs['gene'] in genes):
        found += 1
        if(cds.strand == '+'):
            start = cds.start-args.flank
            stop = cds.stop
        else:
            start = cds.start
            stop = cds.stop + args.flank
        print("\t".join((cds.chrom, str(start), str(stop))))
            
            
sys.stderr.write("\ntotal phage genes: %d\nfound phage genes: %d\n\n" % (found, len(genes)))
        
        
