#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Adds annotation in regard to the antisense transcripts for the already annotated genomic regions (peaks)'''

import argparse
import sys
import os
from collections import defaultdict, Counter
from afbio.sequencetools import coverage2dict


import numpy as np;
import matplotlib.pyplot as plt;
from pybedtools import BedTool

from afbio.pybedtools_af import construct_gff_interval, read_comments


parser = argparse.ArgumentParser(description='Adds annotation in regard to the antisense transcripts for the already annotated genomic regions (peaks)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genomic regions, gff format");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts regions, gff format");
parser.add_argument('--inside', nargs = '?', default=200, type = int, help = "Maximum allowed distance to TSS while inside a gene");
parser.add_argument('--maxd', nargs = '?', default=800, type = int, help = "Maximum allowed distance to TSS");
args = parser.parse_args();


if(os.stat(args.path).st_size == 0):
    sys.exit("###annotate\nInput file is empty, empty output is produced\n")



def annotate_position(peak, tr_dict, maxd, inside):
    if(peak.chrom not in tr_dict):
        return None
    
    center = int(peak.name)
    tr_plus, tr_minus = tr_dict[peak.chrom]
    if(peak.strand == '+'):
        distances = [(tr, center-tr.stop+1) for tr in tr_minus];
    else:
        distances = [(tr, tr.start-center) for tr in tr_plus]
    
    distances = [x for x in distances if x[1]>-1*inside]
    if(not distances):
        return None
    transcript, mindistance = min(distances, key = lambda x: abs(x[1]))

    if(abs(mindistance) <= maxd):
        gtype = 'upstream'
        atg = mindistance + float(transcript.attrs['distance'])
        old_attrs = peak.attrs
        new_attrs = [("annotation", transcript.attrs['annotation']), ("function", transcript.attrs['function']), ("gene", transcript.name), ("genesymbol", transcript.attrs['genesymbol']), ("cg", transcript.attrs.get('cg', 'unknown')), ("tss", mindistance), ("atg", atg), ("gtype", gtype), ("anti", "1")]
        for k, v in new_attrs:
            old_attrs[k] = v;
        

        return construct_gff_interval( peak.chrom, peak.start, peak.stop, 'annotated', score=peak.score, strand=transcript.strand, source='annotate.py', frame='.', attrs=old_attrs.items())
    
    else:
        return None



tr_dict = defaultdict(lambda: ([], []))

for tr in BedTool(args.transcripts):
    if(tr.strand == '+'):
        tr_dict[tr.chrom][0].append(tr)
    else:
        tr_dict[tr.chrom][1].append(tr)

for c in read_comments(args.path):
    print(c)

stat = [0, 0]
for peak in BedTool(args.path):
    peak.attrs['anti'] = '0';
    sys.stdout.write(str(peak))
    anti = annotate_position(peak, tr_dict, args.maxd, args.inside)
    stat[0] += 1
    if(anti):
        sys.stdout.write(str(anti))
        stat[1] += 1;
        

sys.stderr.write("\nTotal number of regions: %d\nNumber of regions with antisense: %d\n\n" % tuple(stat));
    
  
  






























