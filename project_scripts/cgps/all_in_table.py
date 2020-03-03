#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates table with all the information on CgpS time-series project'''

import argparse
import os
import sys

from pybedtools import BedTool, Interval

from afbio.pybedtools_af import read_comments




parser = argparse.ArgumentParser(description='Creates table with all the information on CgpS time-series project');
parser.add_argument('--chap', nargs = '?', required=True, type = str, help = "Path to the chap-seq annotated antisense regions file, regions/all.annotated.gff");
parser.add_argument('--rnaseq', nargs = '?', required=True, type = str, help = "Path to the compiled rnaseq data, gff file");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();



labels = read_comments(args.rnaseq)[0].split("=")[1].split(",")
header = ["Gene ID", "Gene symbol", "Antisense", "Distance ATG", "Distance to TSS"] + ["ChAP T=%s" % x for x in labels] + ["mRNA T=%s" % x for x in labels]


name2rnaseq = dict([ (x.name, x) for x in BedTool(args.rnaseq) ])

res = []
for chap in BedTool(args.chap):
    if(chap.attrs['gtype'] == 'upstream'):
        rnaseq = name2rnaseq[chap.attrs['gene']]
        a = [chap.attrs['gene'], chap.attrs['genesymbol'], chap.attrs['anti'], chap.attrs['atg'], chap.attrs['tss']] + chap.attrs['topcoverage'].split(",") + rnaseq.attrs['expression'].split(":")
        res.append(a);
    
with open(os.path.join(args.outdir, "all_in_table.tsv"), 'w') as f:
    f.write("%s\n" % "\t".join(header))
    for a in res:
        f.write("%s\n" % "\t".join(a))


            
            
