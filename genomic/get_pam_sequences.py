#! /usr/bin/python
'''Finds and reports all PAM sequences for the provided genome'''

import argparse
import sys
import os
from collections import defaultdict

from Bio import SeqIO
from afbio.sequencetools import sliding_window


parser = argparse.ArgumentParser(description='Finds and reports all PAM sequences for the provided genome');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome, fasta format");
#parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format")
args = parser.parse_args();

PAM_LENGTH = 21


def get_pam_sequences(seqrecord, chrname):
    minstart = PAM_LENGTH
    lseq = len(seqrecord)
    
    pam = ("C", "C")

    #forward
    #reverse = seqrecord.reverse_complement()
    #for pos, dn in enumerate(sliding_window(reverse[minstart:], 2), start=minstart):
        #if(dn == pam):
            #start = lseq - pos - 2
            #stop = start + 2 + PAM_LENGTH
            #print(seqrecord[start:stop])
            
    #forward
    for pos, dn in enumerate(sliding_window(seqrecord, 2)):
        if(dn == pam):
            start = pos
            stop = start + 2 + PAM_LENGTH
            seq = str(seqrecord[start:stop])
            if(len(seq) == PAM_LENGTH+2):
                #sys.stderr.write("%s\n" % seq)
                print(">%s|%d|%d|+\n%s" % (chrname, start, stop, seq))
            #print(seqrecord[start:stop])
            
    #reverse
    reverse = seqrecord.reverse_complement()
    for pos, dn in enumerate(sliding_window(reverse, 2)):
        if(dn == pam):
            start = lseq - pos - 2 - PAM_LENGTH
            stop = start + 2 + PAM_LENGTH
            seq = str(seqrecord[start:stop].reverse_complement())
            if(len(seq) == PAM_LENGTH+2):
                print(">%s|%d|%d|-\n%s" % (chrname, start, stop, seq))

    
                
                
for seqrecord in SeqIO.parse(args.path, 'fasta'):
    get_pam_sequences(seqrecord.seq.upper(), seqrecord.name);
