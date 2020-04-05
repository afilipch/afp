#! /usr/local/anaconda3/bin/python
'''Converts custom phage gff annotation to ncbi-like format'''

import argparse
import os
import sys

from pybedtools import BedTool
from afbio.pybedtools_af import construct_gff_interval;

parser = argparse.ArgumentParser(description='Converts ncbi gff annotation to annotated transcripts (custom annotation format)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the ncbi annotation file, gff format");
args = parser.parse_args();





for interval in BedTool(args.path):
    name = 'gene-' + interval.attrs['ID']
    gene = construct_gff_interval(interval.chrom, interval.start, interval.stop, 'gene', score='0', strand=interval.strand, source='custom', frame='.', attrs=[('ID', name), ('gene_biotype', 'protein_coding') ] );
    sys.stdout.write(str(gene))
    cds = construct_gff_interval(interval.chrom, interval.start, interval.stop, 'CDS', score='0', strand=interval.strand, source='custom', frame='.', attrs=[('Parent', name), ('product', interval.attrs.get('product', 'None')), ('Note', interval.attrs.get('inference', 'None')) ] );
    sys.stdout.write(str(cds))
