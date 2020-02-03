#! /usr/bin/python
'''Selects binding peaks according their ability to be counter-silenced by dCas9'''

import argparse
import sys
import os
from collections import defaultdict

from Bio import SeqIO
import numpy as np;
from pybedtools import BedTool
from afbio.pybedtools_af import construct_gff_interval
from afbio.sequencetools import sliding_window


parser = argparse.ArgumentParser(description='Selects binding peaks according their ability to be counter-silenced by dCas9');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks, gff format");
parser.add_argument('--pam', nargs = '?', required=True, type = str, help = "Path to the PAM sequences, fasta format (output of ~/afp/genomic/pam_offset_targets.py)")
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts regions, gff format");
parser.add_argument('--inside', nargs = '?', default=100, type = int, help = "Maximum allowed distance to TSS while inside a gene");
parser.add_argument('--maxd', nargs = '?', default=700, type = int, help = "Maximum allowed distance to TSS");
#parser.add_argument('--pamgap', nargs = '?', default=10, type = int, help = "Minimum allowed gap between PAM motifs");
args = parser.parse_args();

#PAM_LENGTH = 21
HALF_PAM = 10;
PAM_OFFSET = 3


def get_pams_around_peak(center, strand, pam_plus, pam_minus):
    sense = [center - x for x in pam_plus]
    antisense = [x - center for x in pam_minus]
    if(strand == '-'):
        sense, antisense = antisense, sense
        
    return min(sense, key=lambda x: abs(x)), min(antisense, key=lambda x: abs(x))
    



def find_closest(peak, starts_plus, stops_minus, pam_plus, pam_minus, d_threshold, inside):
    center = int(peak.name)
    distances = [('+', x[0], x[1]-center) for x in starts_plus]
    distances.extend([('-', x[0], center-x[1]+1) for x in stops_minus]);
    distances = [x for x in distances if x[2]>-1*inside]
    strand, gene_name, mindistance = min(distances, key = lambda x: abs(x[2]))
    if(abs(mindistance) <= d_threshold):
        pam_sense, pam_antisense = get_pams_around_peak(center, strand, pam_plus, pam_minus)
        pam_min = min([abs(x) for x in [pam_sense, pam_antisense]])
        
        return construct_gff_interval( peak.chrom, peak.start, peak.stop, 'consensus_region', score=str(peak.score), strand=strand, source='un', frame='.', attrs=[("Name", peak.name),("tss", mindistance), ("gene", gene_name), ('pam_sense', pam_sense), ('pam_antisense', pam_antisense), ('pam_min', pam_min) ])
    else:
        return False





###Get PAM centers
pam_plus = [];
pam_minus = [];
for seqrecord in SeqIO.parse(args.pam, 'fasta'):
    chrom, start, stop, strand, off_targets = seqrecord.name.split("|")
    if(off_targets=='None'):
        start = int(start)
        if(strand == '+'):
            pam_plus.append(start+PAM_OFFSET+HALF_PAM) 
        else:
            pam_minus.append(start+HALF_PAM) 
    


tr_list = BedTool(args.transcripts)
starts_plus = [ (x.name, x.start) for x in tr_list if x.strand == '+']
stops_minus = [ (x.name, x.stop) for x in tr_list if x.strand == '-']


for peak in BedTool(args.path):
    annotated = find_closest(peak, starts_plus, stops_minus, pam_plus, pam_minus, args.maxd, args.inside)
    if(annotated):
        sys.stdout.write(str(annotated));

#print([len(x) for x in [seqrecord, pam_minus, pam_plus]])





    






    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
