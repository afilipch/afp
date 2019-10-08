#! /usr/bin/python
'''Finds and reports all PAM sequences for the provided genome'''

import argparse
import sys
import os
from collections import defaultdict
import numpy as np

from Bio import SeqIO
from afbio.sequencetools import sliding_window


parser = argparse.ArgumentParser(description='Analyses the off-targeting of the provided pam sequences');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome, fasta format");
parser.add_argument('--length', nargs = '?', default=12, type = int, help = "Length of the pam sequence to check for the off-targeting");
args = parser.parse_args();

OFFSET = 3
   

def compare_seq(seq, position, first, second, hlen):
    off_targets = [];
    local_sequences = first[seq[:hlen]] + second[seq[hlen:]]
    for cseq, cpos in local_sequences:
        mutations = 0;
        for n1, n2 in zip(seq, cseq):
            mutations += int(n1!=n2);
            if(mutations>1):
                break;
        else:
            if(position!=cpos):
                off_targets.append(cpos)
    return off_targets


sequences = []
for pos, seqrecord in enumerate(SeqIO.parse(args.path, 'fasta')):
    seq = str(seqrecord.seq[OFFSET:OFFSET+args.length])
    sequences.append((seq, pos, seqrecord.name));


hlen = args.length//2    
first = defaultdict(list);
second = defaultdict(list);
for seq, pos, name in sequences:
    first[seq[:hlen]].append((seq, pos));
    second[seq[hlen:]].append((seq, pos));
    

with_offtargets, wo_offtargets = 0, 0
for c, (seq, position, name) in enumerate(sequences, start = 1):
    off_targets = compare_seq(seq, position, first, second, hlen)
    if(off_targets):
        off_targets = ":".join([str(x) for x in off_targets])
        with_offtargets+=1
    else:
        off_targets="None"
        wo_offtargets+=1
    print(">%s|%s\n%s" % (name, off_targets, seq))
    
    if(c % 10000 == 0):
        sys.stderr.write("%d PAM motifs processed\n" % c)
        
        
sys.stderr.write("\nTotal PAM motifs: %d\nWith off-targets: %d\nWithout off-targets: %d" % (with_offtargets + wo_offtargets, with_offtargets, wo_offtargets))


