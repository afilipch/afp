#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Detects AT-rich sequences along the provided genome'''

import argparse
import sys
import os
import numpy as np;
from Bio import SeqIO
from collections import Counter;

from afbio.sequencetools import sliding_window

def get_at_fraction(s):
    return (s.count("A") + s.count('T')) / len(s)

def compare_fraction(seq, minat):
    return get_at_fraction(seq) >= minat

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

    bestvar = max(variants, key = lambda x: calculate_at_extensions(forward, backward, maxgc, x))
    return start-sum([ x[1] for x in backward[:bestvar[1]] ]), start + length + sum([ x[1] for x in forward[:bestvar[0]] ]) 
    
            
        
def get_at_rich_stretches(seq, anchor_length, minlength, maxgc_num, minat_fraction):
    upper_limit = 0
    for position, window in enumerate(sliding_window(seq, anchor_length)):
        if(position >= upper_limit and 'G' not in window and 'C' not in window):
            start, end = get_extensions(seq[upper_limit:], position-upper_limit, anchor_length, maxgc_num, minat_fraction)
            start = max(start, 0)
            if(end-start >= minlength):
                adstart = start + upper_limit
                adend = end + upper_limit
                upper_limit = adend;
                #print(start, end)
                lseq = seq[adstart:adend]
                yield ( adstart, adend, lseq, get_at_fraction(lseq), lseq.count('G') + lseq.count('C') )

                    
def stretch_score(at_fraction, adend, adstart):
    return at_fraction*np.log2(adend - adstart)*10

def stretch_score_flanks(at_fraction, adend, adstart, flanks):
    return at_fraction*np.log2(adend - adstart) + sum(flanks)*0.5




if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(description='Detects AT-rich sequences along the provided genome');
    parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome, fasta format");
    parser.add_argument('--avgat', nargs = '?', required=True, type = float, help = "Average AT content of the provided fasta file (can be found with bin/fastastat.py)");
    parser.add_argument('--anchor', nargs = '?', default=6, type = int, help = "Length of AT only anchor sequence");
    parser.add_argument('--maxgc', nargs = '?', default=2, type = int, help = "Max allowed number of GC inside the stretch");
    parser.add_argument('--minat', nargs = '?', default=0.75, type = float, help = "Min allowed AT content inside the flanks with GC");
    parser.add_argument('--minlength', nargs = '?', default=8, type = int, help = "Min allowed length of a discovered AT stretch");
    parser.add_argument('--flank', nargs = '?', default=0, type = int, help = "Length of the flanks for AT stretches");
    parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
    parser.add_argument('--top', nargs = '?', default=100, type = int, help = "Output top N stretches");
    args = parser.parse_args();
    
    stretches = [];
    for seqrec in SeqIO.parse(args.path, 'fasta'):
        seq = seqrec.seq.upper();
        for adstart, adend, lseq, at_fraction, gc_count in get_at_rich_stretches(seq, args.anchor, args.minlength, args.maxgc, args.minat):
            at_diff = at_fraction - args.avgat
            if(args.flank):
                flanks = seq[max(0, adstart-args.flank) : adstart], seq[adend : adend+args.flank]
                flanks_at = [get_at_fraction(x) - args.avgat for x in flanks]
                score = stretch_score_flanks(at_diff, adend, adstart, flanks_at)
            else:
                flanks_at = 0, 0
                score = stretch_score(at_diff, adend, adstart)
                
            stretches.append(( seqrec.name, adstart, adend, lseq, at_diff*100, gc_count, flanks_at[0]*100, flanks_at[1]*100, score*100))
            
            
    stretches.sort(key = lambda x: x[-1], reverse = True) 
    
    paths = [os.path.join(args.outdir, "at_stretches_top%d_flank%d.%s" % (args.top, args.flank, x)) for x in ('tsv', 'bed')]
    with open(paths[0], 'w') as f1, open(paths[1], 'w') as f2:
        f1.write( "chromosome\tstart\tstop\tseq\tat_diff\tgc_number\tupstream_diff\tdownstream_diff\tscore\n")
        for c, stretch in enumerate(stretches[:args.top], start = 1):
            f2.write( "%s\t%d\t%d\tstretch_%d\t%1.1f\t+\n" % (stretch[0], stretch[1], stretch[2], c, stretch[-1]))
            f1.write( "%s\t%d\t%d\t%s\t%1.1f\t%d\t%1.1f\t%1.1f\t%1.1f\n" % stretch )   



    #testseq2 = "GCGTATATATATATGGAATATGATAAAAAATATATATATAAATAGGAATTAACATATACGTACTCGACTC"
    #print(extend_anchor_forward(testseq2, 30, 6, args.maxgc, args.minat));
    #print(extend_anchor_backward(testseq2, 30, 6, args.maxgc, args.minat));
    #start, end = get_extensions(testseq2, 30, 6, args.maxgc, args.minat)
    #print()
    #print(testseq2)
    #print(testseq2[start:end])
