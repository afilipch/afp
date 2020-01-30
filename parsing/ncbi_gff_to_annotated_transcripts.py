#! /usr/local/anaconda3/bin/python
'''Converts ncbi gff annotation to annotated transcripts (custom annotation format)'''

import argparse
import os
import sys

from pybedtools import BedTool

from afbio.pybedtools_af import construct_gff_interval;

parser = argparse.ArgumentParser(description='Converts ncbi gff annotation to annotated transcripts (custom annotation format)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the ncbi annotation file, gff format");
args = parser.parse_args();


curgene = None;
name2gene = {};
gene_annotation_list = [];

for interval in BedTool(args.path):
    if(interval[2] == 'gene'):
        name2gene[interval.attrs["ID"]] = interval;
    else:
        gene_annotation_list.append(interval)
  
#print(len(name2gene))  

for ga in gene_annotation_list:
    parent_name = ga.attrs.get('Parent')
    if(parent_name):
        parent = name2gene.get(ga.attrs['Parent'])
        if(parent):
            if(parent.strand == '+'):
                distance = ga.start - parent.start
            else:
                distance = parent.stop - ga.stop
            distance = 0;
            aint = construct_gff_interval(parent.chrom, parent.start, parent.stop, 'gene', score='0', strand=parent.strand, source='ncbi_af', frame='.', attrs=[('ID', parent.attrs['Name']), ('genesymbol', parent.attrs.get('gene', 'None')), ('annotation', ga.attrs.get('Note', 'None')), ('product', ga.attrs.get('product', 'None')), ('cds', "%d:%d" % (ga.start+1, ga.stop)), ('tss_variants', '1'), ('distance', str(distance))  ] );
            
            sys.stdout.write(str(aint));

