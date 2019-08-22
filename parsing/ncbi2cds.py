#! /usr/local/anaconda3/bin/python
'''Exctracts coding sequences from ncbi gff annotation file'''

import argparse
import os
import sys

from pybedtools import BedTool

from afbio.pybedtools_af import construct_gff_interval;

parser = argparse.ArgumentParser(description='Exctracts coding sequences from ncbi gff annotation file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the ncbi annotation file, gff format");
args = parser.parse_args();


curgene = None;
cds_list = [];

for interval in BedTool(args.path):
    if(interval[2] == 'gene'):
        curgene = interval;
    if(interval[2] == 'CDS'):
        aint = construct_gff_interval(interval.chrom, interval.start, interval.stop, 'cds', score='0', strand=interval.strand, source='.', frame='.', attrs=[('ID', interval.attrs['Parent']), ('genesymbol', curgene.attrs['Name']), ('annotation', interval.attrs.get('Note', 'None')), ('product', interval.attrs.get('product', 'None'))] );
        cds_list.append(aint)
        
        
cds_list.sort(key=lambda x: x.start);
for aint in cds_list:
    sys.stdout.write(str(aint));
