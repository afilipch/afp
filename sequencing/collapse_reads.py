#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Collapses identical reads'''

import argparse
import sys
import os

from collections import Counter
from math import ceil
from numpy import log10

from afbio.generators import generator_fastq


parser = argparse.ArgumentParser(description='Collapses identical reads for paired-end sequening protocol');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "Path to the sequencing reads");
parser.add_argument('--output', nargs = '?', required=True, type = str, help = "Path to the output files. Only a basename should be provided. For example if \"sample1\" is provided, then output is saved to \"sample1.1.fastq\" and \"sample1.2.fastq\" ");
#parser.add_argument('--ftype', nargs = '?', choices = ['fastq', 'fasta'], default = 'fastq', type = str, help = "Type of the sequences. Can be \'fasta\' or \'fastq\'");
args = parser.parse_args()

#if(args.ftype == 'fastq'):
sequences = [(x[0].seq, x[1].seq) for x in zip(generator_fastq(args.path[0], take = ['seq']), generator_fastq(args.path[1], take = ['seq']))   ]
seq2info = dict([ (x[0], [x[1], [], []])  for x in Counter(sequences).items() ])

def update_quality(qual, current_qual):
    ans = [];
    if(qual):
        for q, c in zip(qual, current_qual):
            if(c>q):
                ans.append(c);
            else:
                ans.append(q);
        ans = "".join(ans);
    else:
        ans = current_qual
    return ans;
        


for read1, read2 in zip(generator_fastq(args.path[0], take = ['seq', 'qual']), generator_fastq(args.path[1], take = ['seq', 'qual'])):
    count, qual1, qual2 = seq2info[read1.seq, read2.seq];
    qual1 = update_quality(qual1, read1.qual);
    qual2 = update_quality(qual2, read2.qual);
    seq2info[read1.seq, read2.seq] = count, qual1, qual2
    


###output the results:
def format_output(seq, qual, count, num):
    return "@read%d_x%d\n%s\n+\n%s\n" % (num+1, count, seq, qual)

with open(args.output + ".1.fastq", 'w') as f1, open(args.output + ".2.fastq", 'w') as f2:    
    for num, (k, v) in enumerate(sorted(seq2info.items(), key = lambda x: x[1][0], reverse=True)):
        seq1, seq2 = k
        count, qual1, qual2 = v;
        f1.write(format_output(seq1, qual1, count, num));
        f2.write(format_output(seq2, qual2, count, num));
        
        
###output additional info:
bins = ['1', '2-10', '11-100', '101-1000', '1001-10000', '>10000']
counts = [x[0] for x in seq2info.values()]
sys.stderr.write("%d reads were collapsed in %d reads\nCompression rate: %1.2f\n\nDistribution of read duplication:\ncollapse factor | number of reads\n" % (sum(counts), len(counts), sum(counts)/len(counts)))
binned = Counter([ceil(log10(x)) for x in counts]);
binned = [x[1] for x in sorted(binned.items(), key=lambda x: x[0])]

for label, num in zip(bins, binned):
    sys.stderr.write("%s:\t%d\n" % (label, num))


