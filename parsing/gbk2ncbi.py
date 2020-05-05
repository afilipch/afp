#! /usr/local/anaconda3/bin/python
'''Converts genebank file into ncbi / MUST BE MADE GENERIC'''

import argparse
import os
import sys
from Bio import SeqIO

from afbio.pybedtools_af import construct_gff_interval;

parser = argparse.ArgumentParser(description='Converts genebank file into gff / MUST BE MADE GENERIC');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genbank file");
args = parser.parse_args();

strand_dict = {1 : '+', -1 : '-'}



#NC_003450.3     RefSeq  gene    4766    5302    .       +       .       ID=gene-NCgl0004;Dbxref=GeneID:1021292;Name=NCgl0004;gbkey=Gene;gene_biotype=protein_coding;gene_synonym=Cgl0005;locus_tag=NCgl0004



for seqrecord in SeqIO.parse(args.path, 'genbank'):
    chrom = seqrecord.name
    chrom = '453-Cg-phage-CL31_S7_L001_R1_001_contig_1'
    for feature in seqrecord.features:
        if(feature.type == 'CDS'):
            temp = feature.qualifiers['label'][0].replace(";", ":")
            a = temp.split(" ")
            geneid = "gene-%s" % a[0]
            product = " ".join(a[1:])
            newint = construct_gff_interval(chrom, feature.location.start, feature.location.end, 'gene', score='0', strand=strand_dict[feature.location.strand], attrs=[('ID', geneid), ('gene_biotype', 'protein_coding'), ('Name', geneid), ('product', product)])
            sys.stdout.write(str(newint))
