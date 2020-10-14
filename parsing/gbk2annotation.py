#! /usr/local/anaconda3/bin/python
'''Converts genebank file into custom annotation gff file'''

import argparse
import os
import sys
from collections import defaultdict

from pybedtools import BedTool
from Bio import SeqIO

from afbio.pybedtools_af import construct_gff_interval;

parser = argparse.ArgumentParser(description='Converts genebank file into custom annotation gff file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genbank file");
args = parser.parse_args();


def print_out(features):
    if(len(features)>2):
        for f in features:
            print(f);
        print()
        print("-"*140)
        print()

strand_dict = {1 : '+', -1 : '-'}
gene_id_dict = defaultdict(list)

for seqrecord in SeqIO.parse(args.path, 'genbank'):
    chrom = seqrecord.name
    cur_features = []
    for feature in seqrecord.features:
        gene_id = feature.qualifiers.get('locus_tag')
        if(gene_id):
            gene_id_dict[(chrom, gene_id[0])].append(feature)

            
for (chrom, geneid), features in gene_id_dict.items():
    fd = dict([ (int(x.type == 'gene'), x) for x in features])
    transcript = fd[0]
    gene = fd[1]
    
    genesymbol = gene.qualifiers.get('gene', [geneid])[0]
    annotation = transcript.qualifiers['product'][0].replace(";", ' ')
    function = transcript.qualifiers.get('note', [''])[0].replace(";", ' ')
    strand = strand_dict[feature.location.strand]
    cds = "%d:%d" % (transcript.location.start, transcript.location.end)
    if(strand == '+'):
        distance = transcript.location.start - gene.location.start
    else:
        distance = gene.location.end - transcript.location.end
    
    newint = construct_gff_interval(chrom, gene.location.start, gene.location.end, transcript.type, score='0', strand=strand, attrs=[('ID', geneid), ('genesymbol', genesymbol), ('annotation', annotation), ('function', function), ('cds', cds), ('distance', distance), ('tss_variants', '1') ])
    sys.stdout.write(str(newint))
    
    
    
    #if('CDS' in [x.type for x in v]):
        
        #print()
        #print("-"*140)
        #print()
        #for el in v:
            #print(el)


