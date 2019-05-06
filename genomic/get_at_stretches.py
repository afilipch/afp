#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Detects AT-rich sequences along the provided genome'''

import argparse
import sys
import os
import numpy as np;
from Bio import SeqIO
from collections import Counter;

from afbio.sequencetools import sliding_window

parser = argparse.ArgumentParser(description='Detects AT-rich sequences along the provided genome');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome, fasta format");
parser.add_argument('--anchor', nargs = '?', default=6, type = int, help = "Length of AT only anchor sequence");
parser.add_argument('--maxgc', nargs = '?', default=2, type = int, help = "Max allowed number of GC inside the stretch");
parser.add_argument('--minat', nargs = '?', default=0.75, type = float, help = "Min allowed AT content inside the flanks with GC");
parser.add_argument('--minlength', nargs = '?', default=8, type = int, help = "Min allowed length of a discovered AT stretch");
args = parser.parse_args();



testseq2 = "GCGTATATATATATGGAATATGATAAAAAATATATATATAAATAGGAATTAACATATACGTACTCGACTC"



def compare_fraction(seq, minat):
    return (seq.count("A") + seq.count('T')) / len(seq) >= minat

def extend_anchor_forward(seq, start, length, maxgc, minat):
    gccount = 0;
    extensions = [];
    current_extension = [];

    for pos, el in enumerate(seq[start+length:]):
        if(el in 'GC'):
            gccount += 1;
            if(gccount>maxgc):
                if(current_extension and compare_fraction(current_extension, minat)):
                    extensions.append(current_extension)
                break;
            else:
                if(current_extension and compare_fraction(current_extension, minat)):
                    extensions.append(current_extension)
                    current_extension = [el];
                else:
                    current_extension.append(el)
        else:
            current_extension.append(el);            
    return extensions;


def extend_anchor_backward(seq, start, length, maxgc, minat):
    gccount = 0;
    extensions = [];
    current_extension = [];

    for pos, el in enumerate(seq[start-1::-1]):
        if(el in 'GC'):
            gccount += 1;
            if(gccount>maxgc):
                if(current_extension and compare_fraction(current_extension, minat)):
                    extensions.append(current_extension)
                break;
            else:
                if(current_extension and compare_fraction(current_extension, minat)):
                    extensions.append(current_extension)
                    current_extension = [el];
                else:
                    current_extension.append(el)
        else:
            current_extension.append(el);            
    return extensions;


def calculate_at_extensions(forward, backward, maxgc, variant):
    total  = forward[:variant[0]] + backward[:variant[1]]
    if(maxgc == sum([x[0] for x in total])):
        nfactor = 0;
    else:
        nfactor = 0.2;
    return 1 - sum([x[0] for x in total])/sum([x[1] for x in total]) - nfactor;


def get_extensions(seq, start, length, maxgc, minat):
    forward = [(x.count("G") + x.count("C"), len(x)) for x in extend_anchor_forward(seq, start, length, maxgc, minat)]
    backward = [(x.count("G") + x.count("C"), len(x)) for x in extend_anchor_backward(seq, start, length, maxgc, minat)]
    #print(forward)
    #print(backward)
    if(not forward and not backward):
        return start, start+length
    elif(not forward):
        return start-sum([x[1] for x in backward]), start + length
    elif(not backward):
        return start, start + length + sum([ x[1] for x in forward ])
    variants = [];
    for fc in range(len(forward)+1):
        gclimit = maxgc - sum([x[0] for x in forward[:fc]])
        #print()
        #print(gclimit)
        #print()
        for bc in range(len(backward)+1):
            curgc = sum([x[0] for x in backward[:bc]])
            #print(curgc);
            if(curgc>gclimit):
                variants.append((fc, bc-1));
                break;
        else:
            variants.append((fc, len(backward)))
            
    #print(variants)
    #print([calculate_at_extensions(forward, backward, maxgc, x) for x in variants])
    bestvar = max(variants, key = lambda x: calculate_at_extensions(forward, backward, maxgc, x))
    #print(bestvar)
    return start-sum([ x[1] for x in backward[:bestvar[1]] ]), start + length + sum([ x[1] for x in forward[:bestvar[0]] ]) 
    
            
        
    
    
    

#print(extend_anchor_forward(testseq2, 30, 6, args.maxgc, args.minat));
#print(extend_anchor_backward(testseq2, 30, 6, args.maxgc, args.minat));
#start, end = get_extensions(testseq2, 30, 6, args.maxgc, args.minat)
#print()
#print(testseq2)
#print(testseq2[start:end])
                    
        


upper_limit = 0;
for seqrec in SeqIO.parse(args.path, 'fasta'):
    seq = seqrec.seq.upper();
    for position, window in enumerate(sliding_window(seq, args.anchor)):
        if(position >= upper_limit and 'G' not in window and 'C' not in window):
            start, end = get_extensions(seq[upper_limit:], position-upper_limit, args.anchor, args.maxgc, args.minat)
            if(end-start >= args.minlength):
                adstart = start + upper_limit
                adend = end + upper_limit
                upper_limit = adend;
                lseq = seq[adstart:adend]
                print( "%s\t%d\t%d\tr%d_%d\t%d\t+" % ('chr1', adstart, adend, adstart, adend, lseq.count('G') + lseq.count('C')))









        
