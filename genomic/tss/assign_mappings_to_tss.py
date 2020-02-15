#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Assigns mappings to transcripts with multiple TSS'''

import argparse
import os
import sys
from collections import defaultdict

import numpy as np;
import matplotlib.pyplot as plt;

from pybedtools import BedTool

from afbio.sequencetools import coverage2dict
from afbio.generators import get_only_files


parser = argparse.ArgumentParser(description='Assigns mappings to transcripts with multiple TSS');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the coverage folder, bed format");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts with multiple tss (custom made), gff format");
parser.add_argument('--order', nargs = '+', required = True, type = str, help = "Labels are output in this order");
parser.add_argument('--outdir', nargs = '?', required = True, type = str, help = "Path to the output directory");
args = parser.parse_args();

def convert_to_tpm(cds2cov):
    scaling_factor = sum([x for x in cds2cov.values()])/(10**6)
    res = dict([ (x[0], x[1]/scaling_factor) for x in cds2cov.items() ])
    return res;


def find_fraction(e1, l1, e2, l2):
    return ((e2-e1)/(l2-l1))/(e1/l1)
    

def tss_to_fraction(expr_list, variants):
    expr_list.sort()
    l_list = list(sorted( [x[2] - x[1] for x in variants] ))
    if(variants[0][0] == '+'):
        tss_list = list(sorted( [x[1] for x in variants], reverse=True));
    else:
        tss_list = list(sorted( [x[2] for x in variants]));
        
    fractions = [];
    curfraction = 1;
    for (e1, l1), (e2, l2) in zip( zip(expr_list, l_list), zip(expr_list[1:], l_list[1:]) ):
        if(not e1):
            fractions = [0]*len(l_list)
            break;
        fraction = find_fraction(e1, l1, e2, l2)
        #print(fraction)
        if(fraction<=1):
            fractions.append(curfraction*(1-fraction))
            curfraction *= fraction
        else:
            fractions.append(1)
            curfraction *= fraction
    else:
        fractions.append(curfraction)
            
     
    return fractions, tss_list




files = get_only_files(args.path);
plus_dict = defaultdict(list)
for pf in sorted([x for x in files if "plus" in x]):
    plus_dict[os.path.basename(pf).split("_")[0]].append(pf)
minus_dict = defaultdict(list)
for pf in sorted([x for x in files if "minus" in x]):
    minus_dict[os.path.basename(pf).split("_")[0]].append(pf)

cov_list = [];
for label in args.order:
    cov_list.append( list(zip(plus_dict[label], minus_dict[label])) )


transcripts = BedTool(args.transcripts)
sample_res_list = defaultdict(lambda: [[] for _ in args.order])
for sample, cl in enumerate(cov_list):
    for plus_cov, minus_cov in cl:
        print(plus_cov, minus_cov)
        covdict = {'+': coverage2dict(plus_cov), '-': coverage2dict(minus_cov)}
        tss2cov = defaultdict(list)
        cds2variants = defaultdict(list)
        cds2cov = {}
        for transcript in transcripts:
            chrom_cov = covdict[transcript.strand][transcript.chrom]
            total_coverage =  sum(chrom_cov[transcript.start:transcript.stop])
            tss2cov[transcript.name].append(total_coverage)
            cds2variants[transcript.name].append((transcript.strand, transcript.start, transcript.stop))
            
            start, stop = [int(x) for x in transcript.attrs['cds'].split(":")]
            rbp = np.mean(chrom_cov[start:stop])
            cds2cov[transcript.name] = rbp
            
        cds2cov = convert_to_tpm(cds2cov);
            
        for tname, expr_list in tss2cov.items():
            variants = cds2variants[tname]
            if(len(variants)>1):
                fractions, tss_list = tss_to_fraction(expr_list, variants)
                sample_res_list[tname][sample].append((  fractions, tss_list, cds2cov[tname] ))
                
        #for tname, el in sample_res_list.items():
            #print(tname, el);
        
   
header = "sample\t" + "\t".join(args.order)
for tname, el in sample_res_list.items():
    with open(os.path.join(args.outdir, "%s.tsv" % tname), 'w') as f:
        #print(el)
        f.write("%s\n" % header)
        f.write("expression\t%s\n" % "\t".join([",".join(["%1.2f" % y[2] for y in x]) for x in el]))
        for vpos, variant in enumerate(el[0][0][1]):
            #print(variant)
            #print( "expression\t%s\n" % "\t".join([",".join(["%1.4f" % y[0][vpos] for y in x]) for x in el]) )
            f.write("%s\t%s\n" % (variant, "\t".join([",".join(["%1.4f" % y[0][vpos] for y in x]) for x in el]))  )
        

        
        
        
#scaling_factor = sum([x for x in transcript2cov.values()])/(10**6)
#for transcript in transcripts:
    #if(transcript[2] in ['gene', 'pseudogene'] and transcript.attrs['gene_biotype'] == 'protein_coding'): 
        #transcript.attrs['tpm'] = '%1.3f' % (transcript2cov[transcript.name]/scaling_factor)
        #sys.stdout.write(str(transcript))
        
        























