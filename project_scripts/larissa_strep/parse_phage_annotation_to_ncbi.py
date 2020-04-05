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

#NC_003888.3     RefSeq  gene    14044   14454   .       -       .       ID=gene-SCO0012;Dbxref=GeneID:1095452;Name=SCO0012;gbkey=Gene;gene_biotype=protein_coding;gene_synonym=SCJ30.07c;locus_tag=SCO0012

#NC_003888.3     RefSeq  CDS     14044   14454   .       -       0       ID=cds-NP_624373.1;Parent=gene-SCO0012;Dbxref=Genbank:NP_624373.1,GeneID:1095452;Name=NP_624373.1;Note=SCJ30.07c%2C unknown%2C len: 136 aa%3B predicted by GC Frameplot%2C Hidden Markov model and amino acid usage.;gbkey=CDS;locus_tag=SCO0012;product=hypothetical protein;protein_id=NP_624373.1;transl_table=11



# phiSco2_Bielefeld       Prodigal:2.6    CDS     1400    2647    .       -       0       ID=phiSco2_00003;Parent=phiSco2_00003_gene;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:Viruses2.aa:YP_009208303.1;locus_tag=phiSco2_00003;product=hypothetical protein SEA_AMELA_25 



for interval in BedTool(args.path):
    if(interval[2] == 'CDS'):
        name = 'gene-' + interval.attrs['ID']
        gene = construct_gff_interval(interval.chrom, interval.start, interval.stop, 'gene', score='0', strand=interval.strand, source='custom', frame='.', attrs=[('ID', name), ('gene_biotype', 'protein_coding') ] );
        sys.stdout.write(str(gene))
        cds = construct_gff_interval(interval.chrom, interval.start, interval.stop, 'CDS', score='0', strand=interval.strand, source='custom', frame='.', attrs=[('Parent', name), ('product', interval.attrs.get('product', 'None')), ('Note', interval.attrs.get('inference', 'None')) ] );
        sys.stdout.write(str(cds))
                                 
    
