#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Converts annotated peaks gff file into tsv file'''

import argparse
import sys, os

from pybedtools import BedTool
from afbio.pybedtools_af import read_comments
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Converts annotated peaks gff file into tsv file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated peaks, gff");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file")
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to output directory")
parser.add_argument('--flanks', nargs = '+', default=[30], type = int, help = "Sizes of the flanks for the peak sequences")
args = parser.parse_args();


def interval2seq(interval, reference, flank):
    start = int(interval.name) - flank
    stop = int(interval.name) + flank +1
    if(interval.strand == '+'):
        return str(reference[interval.chrom][start:stop].seq.upper())
    elif(interval.strand == '-'):
        return str(reference[interval.chrom][start:stop].seq.reverse_complement().upper())


labels  = read_comments(args.path)[0][1:].strip().split(",")
genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))  
regions = BedTool(args.path)
  
flank_header = ["Seq (%d bp)" % (x*2+1) for x in args.flanks]
attrs = ['gene', 'genesymbol', 'annotation', 'function', 'tss', 'atg', 'gtype', 'anti_gene', 'anti_genesymbol', 'anti_tss']
if('cg' in regions[0].attrs.keys()):
    attrs.insert(2, 'cg')
header = ['Chrom', 'Start', 'End', 'Strand'] + flank_header + labels + attrs


with open(os.path.join(args.outdir, "annotated.tsv"), 'w') as f1, open(os.path.join(args.outdir, "annotated.fa"), 'w') as f2:
    f1.write("\t".join(header) + "\n")
    for interval in regions:
        seqs = [interval2seq(interval, genome, x) for x in args.flanks]
        f1.write( "\t".join( [str(x) for x in (interval.chrom, interval.start, interval.end, interval.strand)] + seqs + interval.attrs["topcoverage"].split(',') + [interval.attrs[x] for x in attrs]) + "\n")
        f2.write( ">%s(%s):%d:%d (%s)\n%s\n" % (interval.chrom, interval.strand, interval.start, interval.end, interval.attrs['gene'], seqs[-1]))
