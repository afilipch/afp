#! /usr/bin/python
'''Adds an information to transcripts whether they belong to phage or not'''

import argparse
import sys
import os

from pybedtools import BedTool;
from afbio.pybedtools_af import intersection2gff;

parser = argparse.ArgumentParser(description='Adds an information to transcripts whether they belong to phage or not');
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the annotated transcripts, gff file");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
args = parser.parse_args();


phages = BedTool(args.phages);
transcripts = BedTool(args.transcripts);
for inter in transcripts.intersect(b=phages, wao =True, f = 0.5):
    i1, i2 = intersection2gff(inter);
    if(i2):
        i1.attrs['phage'] = '1'
    else:
        i1.attrs['phage'] = '0'
    sys.stdout.write(str(i1));
        
    
