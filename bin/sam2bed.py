#! /usr/local/anaconda3/bin/python
'''Converts sam records with multimappers into bedfile'''

import argparse
import os
import sys
import numpy as np;
import matplotlib.pyplot as plt;

import pysam
from Bio.Seq import reverse_complement

from afbio.generators import generator_doublesam

parser = argparse.ArgumentParser(description='Converts sam records with multimappers into bedfile');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the sam file");
parser.add_argument('--multimappers', nargs = '?', default='', type = str, help = "Path to store multimapped reads. If not set, multimapped reads are discrarded");
args = parser.parse_args();


def tobed(s1, s2, count=0, norm=0):
    score = s1.get_tag('AS') + s2.get_tag('AS')
    if(s1.is_reverse):
        strand = '+'
    else:
        strand = '-'
    start = min(s1.reference_start, s2.reference_start);
    end = max(s1.reference_end, s2.reference_end);
    chrom = s1.reference_name
    if(count):
        name = '_'.join((s1.query_name, str(count), str(norm)))
    else:
        name = s1.query_name
    return '\t'.join( [str(x) for x in (chrom, start, end, name, score, strand)] )
    


def convert_mappings_multi(readmappings):
    norm = len(readmappings);
    ans = [];
    for count, (s1, s2) in enumerate(readmappings):
        ans.append(tobed(s1, s2, count+1, norm))
    return ans;

def convert_mappings(readmappings):
    return tobed(*readmappings[0]);
            
if(args.multimappers):
    fmult = open(args.multimappers, 'w');



samfile = pysam.AlignmentFile(args.path)
curname = '';
readmappings = [];
for seg1, seg2 in generator_doublesam(samfile):
    if(seg1.query_name != seg2.query_name):
        sys.stderr.write("WARNING: Consequtive paired reads have different names:\t%s\t%s\n" % (seg1.query_name, seg2.query_name))
    if(seg1.query_name == curname):
        readmappings.append((seg1, seg2))
    else:
        if(readmappings):
            if(len(readmappings)==1):
                print(convert_mappings(readmappings))
            elif(args.multimappers):
                for el in convert_mappings_multi(readmappings):
                    fmult.write(el+"\n")
        readmappings = [(seg1, seg2)]
        curname = seg1.query_name

else:
    if(len(readmappings)==1):
        print(convert_mappings(readmappings))
    elif(args.multimappers):
        for el in convert_mappings_multi(readmappings):
            fmult.write(el+"\n")
            
if(args.multimappers):
    fmult.close();            
        
