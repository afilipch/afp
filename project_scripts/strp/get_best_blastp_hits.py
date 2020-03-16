#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Selects the best hits from blastp output'''

import argparse
import os
import sys

from Bio import SearchIO


parser = argparse.ArgumentParser(description='Selects the best hits from blastp output');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to blast output");
#parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Path to the table with proteins from human genome atlas");
parser.add_argument('--evalue', nargs = '?', default=10**-20, type = float, help = "minimum allowed e-value");
args = parser.parse_args()


#def process_blast_hit(curhit):
    #if(len(curhit)>2):
        #print(curhit);
        #print()
        #print()

#def read_blast(path):
    #hits = [];
    #curhit = [];
    #with open(path) as f:
        #for l in f:
            #if(l.startswith("Query=")):
                #if(curhit):
                    #hit = process_blast_hit(curhit)
                    #if(hit):
                        #hits.append(hit)
                #curhit = [l];
            #else:
                #curhit.append(l)
                
            

#read_blast(args.path)
#custom_fields = ['qseqid', 'sseqid', 'evalue', 'score']
#custom_fields = 'qseqid sseqid evalue score'
for qresult in SearchIO.parse(args.path, 'blast-tab'):
    #print(len(qresult))
    for hit in qresult:
        print(len(hit))
        #if(len(hit)>1):
            #for local_hit in hit:
                #print(local_hit)
                #print()
            #print("-"*120)
