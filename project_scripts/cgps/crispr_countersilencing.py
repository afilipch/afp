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



def find_closest(peak, starts_plus, stops_minus, d_threshold, inside):
    center = int(peak.name)
    distances = [('+', x[0], x[1]-center) for x in starts_plus]
    distances.extend([('-', x[0], center-x[1]+1) for x in stops_minus]);
    distances = [x for x in distances if x[2]>-1*inside]
    strand, gene_name, mindistance = min(distances, key = lambda x: abs(x[2]))
    if(abs(mindistance) <= d_threshold):
        return construct_gff_interval( peak.chrom, peak.start, peak.stop, 'consensus_region', score='0', strand=strand, source='un', frame='.', attrs=[("tss", mindistance)] )


def get_pam_sequences(seqrecord, mingap):
    minstart = 21
    pam = ("G", "G")
    res = [];
    for genseq in (seqrecord, seqrecord.reverse_complement()):
        pamlist = [0]*21
        curgap = mingap;
        
        for pos, dn in enumerate(sliding_window(genseq[minstart:], 2), start=minstart):
            if(dn == pam and curgap>=mingap):
                pamlist.append(1)
                curgap = 0;
            else:
                pamlist.append(0)
                curgap+=1;
        res.append(pamlist);
    return res
    



#tr_list = BedTool(args.transcripts)
#starts_plus = [ (x.name, x.start) for x in tr_list if x.strand == '+']
#stops_minus = [ (x.name, x.stop) for x in tr_list if x.strand == '-']

seqrecord = next(SeqIO.parse(args.genome, 'fasta')).seq.upper()[:200]
pam_plus, pam_minus = get_pam_sequences(seqrecord, args.pamgap);





    






    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
