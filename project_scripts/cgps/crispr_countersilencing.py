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
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format")
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts regions, gff format");
parser.add_argument('--inside', nargs = '?', default=100, type = int, help = "Maximum allowed distance to TSS while inside a gene");
parser.add_argument('--maxd', nargs = '?', default=700, type = int, help = "Maximum allowed distance to TSS");
parser.add_argument('--pamgap', nargs = '?', default=10, type = int, help = "Minimum allowed gap between PAM motifs");
args = parser.parse_args();

PAM_LENGTH = 21
HALF_PAM = 10;


def get_pams_around_peak(center, strand, pam_plus, pam_minus, flank, length):
    start = max(center - flank, 0)
    stop = center + flank
    local_plus = pam_plus[start:stop];
    local_minus = pam_minus[start:stop];


    sense = [x[0] + 2 + length - flank for x in enumerate(local_plus) if x[1]==1]
    antisense = [flank - x[0] - 2 + length for x in enumerate(local_minus) if x[1]==1]
    sense.append(flank*2)
    antisense.append(flank*2)
    if(strand == '-'):
        sense, antisense = antisense, sense
        
    return min(sense, key=lambda x: abs(x)), min(antisense, key=lambda x: abs(x))
    



def find_closest(peak, starts_plus, stops_minus, pam_plus, pam_minus, d_threshold, inside, pam_flank, pam_half_length):
    center = int(peak.name)
    distances = [('+', x[0], x[1]-center) for x in starts_plus]
    distances.extend([('-', x[0], center-x[1]+1) for x in stops_minus]);
    distances = [x for x in distances if x[2]>-1*inside]
    strand, gene_name, mindistance = min(distances, key = lambda x: abs(x[2]))
    if(abs(mindistance) <= d_threshold):
        pam_sense, pam_antisense = get_pams_around_peak(center, strand, pam_plus, pam_minus, pam_flank, pam_half_length)
        pam_min = min([abs(x) for x in [pam_sense, pam_antisense]])
        
        return construct_gff_interval( peak.chrom, peak.start, peak.stop, 'consensus_region', score=str(peak.score), strand=strand, source='un', frame='.', attrs=[("Name", peak.name),("tss", mindistance), ("gene", gene_name), ('pam_sense', pam_sense), ('pam_antisense', pam_antisense), ('pam_min', pam_min) ])
    else:
        return False


def get_pam_sequences(seqrecord, mingap, minstart):
    pam = ("G", "G")
    res = [];
    for genseq in (seqrecord.reverse_complement(), seqrecord):
        pamlist = [0]*(minstart+1)
        curgap = mingap;
        
        for pos, dn in enumerate(sliding_window(genseq[minstart:], 2), start=minstart):
            if(dn == pam and curgap>=mingap):
                pamlist.append(1)
                curgap = 0;
            else:
                pamlist.append(0)
                curgap+=1;
        res.append(pamlist);
    res[0] = res[0][::-1]
    return res
    



tr_list = BedTool(args.transcripts)
starts_plus = [ (x.name, x.start) for x in tr_list if x.strand == '+']
stops_minus = [ (x.name, x.stop) for x in tr_list if x.strand == '-']

seqrecord = next(SeqIO.parse(args.genome, 'fasta')).seq.upper()
pam_plus, pam_minus = get_pam_sequences(seqrecord, args.pamgap, PAM_LENGTH);

for peak in BedTool(args.path):
    annotated = find_closest(peak, starts_plus, stops_minus, pam_plus, pam_minus, args.maxd, args.inside, 300, HALF_PAM)
    if(annotated):
        sys.stdout.write(str(annotated));

#print([len(x) for x in [seqrecord, pam_minus, pam_plus]])





    






    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
