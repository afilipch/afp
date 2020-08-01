#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Collapses identical reads'''

import argparse
import sys
import os

from collections import Counter
from math import ceil
from numpy import log10

from afbio.generators import generator_fastq


parser = argparse.ArgumentParser(description='Collapses identical reads for single-end sequening protocol');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the sequencing reads");
parser.add_argument('--output', nargs = '?', required=True, type = str, help = "Path to the output files. Only a basename should be provided. For example if \"sample1\" is provided, then output is saved to \"sample1.fastq\" ");
parser.add_argument('--log', nargs = '?', type = str, help = "Path to the log file");
args = parser.parse_args()

#if(args.ftype == 'fastq'):
sequences = [x.seq for x in generator_fastq(args.path, take = ['seq']) ]
seq2info = dict([ (x[0], [x[1], []])  for x in Counter(sequences).items() ])

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
        


for read in generator_fastq(args.path, take = ['seq', 'qual']):
    count, qual = seq2info[read.seq];
    qual = update_quality(qual, read.qual);
    seq2info[read.seq] = count, qual
    


###output the results:
def format_output(seq, qual, count, num):
    return "@read%d_x%d\n%s\n+\n%s\n" % (num+1, count, seq, qual)

with open(args.output + ".fastq", 'w') as f:    
    for num, (k, v) in enumerate(sorted(seq2info.items(), key = lambda x: x[1][0], reverse=True)):
        seq = k
        count, qual = v;
        f.write(format_output(seq, qual, count, num));

        
###output additional info:
bins = ['1', '2-10', '11-100', '101-1000', '1001-10000', '>10000']
counts = [x[0] for x in seq2info.values()]
binned = Counter([ceil(log10(x)) for x in counts]);
binned = [x[1] for x in sorted(binned.items(), key=lambda x: x[0])]

stat = ["%d reads were collapsed in %d reads\nCompression rate: %1.2f\n\nDistribution of read duplication:\ncollapse factor | number of reads\n" % (sum(counts), len(counts), sum(counts)/len(counts))]
for label, num in zip(bins, binned):
    stat.append("%s:\t%d\n" % (label, num))
    
if(args.log):
    with open(args.log, 'w') as f:
        f.write("\nFile %s is processed\n\n" % args.path)
        f.write("\n".join(stat) + "\n" + "-"*120)
sys.stderr.write("\n".join(stat) + "\n\n" + "-"*120)


