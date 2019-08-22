#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Annotates the table of differentially expressed genes'''

import argparse
import os
import sys

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Annotates the table of differentially expressed genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the differentially expressed genes, tsv table");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the gene annotation, gff format");
args = parser.parse_args();

geneid2annotation = {};
for interval in BedTool(args.annotation):
    geneid2annotation[interval.attrs['geneid']] = (interval.name, interval.attrs['annotation'], interval.attrs['function'])
    
with open(args.path) as f:
    header = next(f).strip().split("\t")
    header[1] = header[1] + "(TPM)"
    header[2] = header[2] + "(TPM)"
    header = ["name"] + header + ["gene_annotation", "gene_function"]
    print("\t".join(header))
    for l in f:
        a = l.strip().split("\t")
        name, annotation, function = geneid2annotation.get(a[0], (a[0], 'unknown', 'unknown'));
        a = [name] + a + [annotation, function]
        print("\t".join(a))
        
        

