#! /usr/local/anaconda3/bin/python
'''Exctracts coding sequences from ncbi gff annotation file'''

import argparse
import os
import sys

from pybedtools import BedTool
from Bio import SeqIO

from afbio.pybedtools_af import construct_gff_interval;

parser = argparse.ArgumentParser(description='Exctracts coding sequences from ncbi gff annotation file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the ncbi annotation file, gff format");
args = parser.parse_args();

strand_dict = {1 : '+', -1 : '-'}

for seqrecord in SeqIO.parse(args.path, 'genbank'):
    chrom = seqrecord.name
    for feature in seqrecord.features:
        if(feature.type == 'gene'):
            geneid = feature.qualifiers['locus_tag']
            genesymbol = feature.qualifiers.get('gene', geneid)[0]
            geneid = geneid[0]
            print("%s\t%d\t%d\t%s\t0\t%s" % (chrom, feature.location.start, feature.location.end, genesymbol, strand_dict[feature.location.strand]) )


