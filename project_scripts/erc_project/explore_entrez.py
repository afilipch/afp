#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explore prophage database'''

import argparse
import os
import sys
#import copy
from collections import defaultdict, Counter
from itertools import product;
from glob import glob


from Bio import Entrez
import numpy as np;
#from scipy.stats import pearsonr
#import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Explore prophage database');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "UNIX-like pattern path to the prophage features files");
parser.add_argument('--failed', nargs = '?', required=True, type = str, help = "Path to the recovery file. Accessions which failed to be processed will be stored there");
parser.add_argument('--add_failed', nargs = '?', type = str, help = "Path to the accession numbers which failed to be downloaded. If set only these numbers will be processed");
#parser.add_argument('--maxshift', nargs = '?', default=50, type = int, help = "Max allowed shift (in nucleotides) of the peak top position downstream to start of the gene, to be still counted as peak upstream the gene");
args = parser.parse_args();

Entrez.email = "vassenego@mail.ru"


if(args.add_failed):
    add_failed = set()
    with open(args.add_failed) as f:
        for l in f:
            add_failed.add(l.strip());
    for bacteria_accession in add_failed:
        handle = Entrez.efetch(db='nucleotide', id = bacteria_accession, rettype = 'gb', retmode = 'xml')
        record = Entrez.read(handle)[0]
        print("%s\t%s\t%s" % (bacteria_accession, record['GBSeq_taxonomy'], record['GBSeq_length']));
        
            
else:
    already_processed_bacteria = set();            
    with open(args.failed, 'w') as failed:
        for c, fpath in enumerate(glob(args.path)[:]):
            bacteria_accession = os.path.basename(fpath).split("_")[0]
            if(bacteria_accession not in already_processed_bacteria):
                try:
                    handle = Entrez.efetch(db='nucleotide', id = bacteria_accession, rettype = 'gb', retmode = 'xml')
                    record = Entrez.read(handle)[0]
                    print("%s\t%s\t%s" % (bacteria_accession, record['GBSeq_taxonomy'], record['GBSeq_length']));
                    already_processed_bacteria.add(bacteria_accession)
                except:
                    sys.stderr.write("\n\nFailed to process %s\n\n" % bacteria_accession);
                    failed.write("%s\n" % bacteria_accession)
                
            if(c and c % 100 == 0):
                sys.stderr.write("%d files processed\n" % c);
