#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Generates genomic coverage based on provided mappings (paired sam/bam files). This script also outputs miscelaneous statistics for the mapped paired-end reads'''
import argparse
import os
import sys
#sys.stderr.write("%s\n\n" % sys.version)
import numpy as np;
import matplotlib.pyplot as plt;
from collections import defaultdict
from itertools import product

import pysam
from Bio.Seq import reverse_complement
from Bio import SeqIO
from matplotlib import pyplot as plt;

from afbio.generators import generator_single_mappings, generator_paired_mappings

parser = argparse.ArgumentParser(description='Converts sam records with multimappers into bedfile');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the sam file");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the reference genome, fasta format");
parser.add_argument('--outstat', nargs = '?', required=True, type = str, help = "Path to a folder for statistics");
parser.add_argument('--outcoverage', nargs = '?', required=True, type = str, help = "Path to a folder for coverage");
parser.add_argument('--ambiguous', nargs = '?', default=0, const=1, type = int, help = "If set ambiguous mappings will be also counted");
parser.add_argument('--paired', nargs = '?', default=False, const=True, type = bool, help = "has to be set if the mappings arise from paired-end reads");
parser.add_argument('--collapsed', nargs = '?', default=False, const=True, type = bool, help = "has to be set to count for collapsed reads");
#parser.add_argument('--multimappers', nargs = '?', default='', type = str, help = "Path to store multimapped reads. If not set, multimapped reads are discrarded");
args = parser.parse_args();

def get_stat_paired(pair):

    if(pair[0].is_reverse):
        strand = '+'
    else:
        strand = '-'
    
    chrom = pair[0].reference_name
    start = min([x.reference_start for x in pair]);
    end = max([x.reference_end for x in pair]);
    score = sum([int(x.get_tag('AS')) for x in pair])
            
    return chrom, start, end, score, strand 


def get_stat_single(s):
    if(s.is_reverse):
        strand = '+'
    else:
        strand = '-'
    
    chrom = s.reference_name
    start = s.reference_start 
    end = s.reference_end
    score = int(s.get_tag('AS'))
            
    return chrom, start, end, score, strand 
    
    
def readmappings2stat(readmappings, func, ambiguous):
    if(not readmappings):
        return []
    if(not ambiguous and len(readmappings)>1):
            return []
    return [func(x) for x in readmappings]


if(args.paired):
    generator_mappings = generator_paired_mappings;
    get_stat = get_stat_paired
else:
    generator_mappings = generator_single_mappings;
    get_stat = get_stat_single



scores = defaultdict(int)
fragment_lengthes = defaultdict(int)
coverage = {};
starts = {};
ends = {};
for seqrecord in SeqIO.parse(args.genome, 'fasta'):
    coverage[(seqrecord.name, '+')] = np.zeros(len(seqrecord));
    coverage[(seqrecord.name, '-')] = np.zeros(len(seqrecord));
    starts[(seqrecord.name, '+')] = np.zeros(len(seqrecord));
    starts[(seqrecord.name, '-')] = np.zeros(len(seqrecord));
    ends[(seqrecord.name, '+')] = np.zeros(len(seqrecord)+1);
    ends[(seqrecord.name, '-')] = np.zeros(len(seqrecord)+1);


samfile = pysam.AlignmentFile(args.path)
for readmappings in generator_mappings(samfile):
    stat_list = readmappings2stat(readmappings, get_stat, args.ambiguous)
    if(stat_list):
        count = 1.0/len(stat_list)
        if(args.collapsed):
            count *= int(readmappings[0].query_name.split("_")[-1][1:] )
        for chrom, start, end, score, strand in stat_list:
            coverage[chrom, strand][start:end] += count;
            fragment_lengthes[end-start] += count;
            starts[chrom, strand][start] += count;
            ends[chrom, strand][end] += count;
            scores[score] += count;
            
    


### GET SAMPLE basename
basename = ".".join((os.path.basename(args.path)).split(".")[:-1]);





#####################################################################################################################################
###RNAse cut preferences
s_nucls = defaultdict(int)
e_nucls = defaultdict(int)
#print(starts)
for seqrecord in SeqIO.parse(args.genome, 'fasta'):
    s_plus = starts[(seqrecord.name, '+')];
    #print('bu')
    s_minus = starts[(seqrecord.name, '-')];
    e_plus = ends[(seqrecord.name, '+')];
    e_minus = ends[(seqrecord.name, '-')];
    for s, e, n1, n2 in zip(s_plus[1:], e_plus[1:], seqrecord.seq, seqrecord.seq[1:]):
        nn = n1+n2
        s_nucls[nn] += s;
        e_nucls[nn] += e;
        
    for s, e, n1, n2 in zip(s_minus[1:], e_minus[1:], seqrecord.seq, seqrecord.seq[1:]):
        nn = reverse_complement(n1+n2)
        s_nucls[nn] += e;
        e_nucls[nn] += s;

def rnase_plot(ndict, output, title, normed=False):
    xticklabels = ["".join(x) for x in product('ACTG', repeat=2)]
    bars = np.array([ndict[x] for x in xticklabels])
    brange = range(16)
    if(normed):
        bars = bars/sum(bars);
        ylabel = 'Fraction'
    else:
        ylabel = 'Counts'
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.bar(brange, bars, 0.5, color='lightblue')
    plt.title(title)
    plt.xticks(brange, ["|".join(x) for x in xticklabels])
    plt.ylabel(ylabel)
    plt.savefig(output, format='png');
    plt.clf()
    return xticklabels, bars
    
rnase_plot(s_nucls, os.path.join(args.outstat, '%s.rnase.5prime.png' % basename), "%s 5'end cut frequences" % basename, normed=True)
rnase_plot(e_nucls, os.path.join(args.outstat, '%s.rnase.3prime.png' % basename), "%s 3'end cut frequences" % basename, normed=True)
        



#####################################################################################################################################
###Length distribution
def draw_small_scale_distr(distr, output, title, xlabel, normed = False, fontsize=24):
    feature_vals = distr.keys();
    lower = min(feature_vals)
    upper = max(feature_vals)+1
    brange = range(lower, upper)
    bars = np.array([distr.get(x,0) for x in brange])
    if(normed):
        bars = bars/sum(bars);
        ylabel = 'Fraction'
    else:
        ylabel = 'Counts'
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.bar(brange, bars, 1, color='lightblue')
    #plt.title(title)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(output, format='png');
    plt.clf()
    return list(brange), bars
    
draw_small_scale_distr(fragment_lengthes, os.path.join(args.outstat, '%s.fragment_lengthes.png' % basename), "%s fragment length distribution" % basename, 'length', normed=True)
draw_small_scale_distr(scores, os.path.join(args.outstat, '%s.scores.png' % basename), "%s alignment score distribution" % basename, 'score', normed=True)


#####################################################################################################################################
###Genomic coverage
with open(os.path.join(args.outcoverage, '%s.plus.bed' % basename), 'w') as fplus, open(os.path.join(args.outcoverage, '%s.minus.bed' % basename), 'w') as fminus:
    for (chrom, strand), localcov in coverage.items():
        sys.stderr.write("Total coverage of chromosome %s on strand %s is %d\n"  % (chrom, strand, sum(localcov)))
        if(strand == '+'):
            for pos, cov in enumerate(localcov): 
                fplus.write("%s\t%d\t%d\n" % (chrom, pos, cov));
        else:
            for pos, cov in enumerate(localcov): 
                fminus.write("%s\t%d\t%d\n" % (chrom, pos, cov));
                

    
        
    
        
        
        
        
        
        
        
        
       
        
