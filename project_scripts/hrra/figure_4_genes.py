#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Shows an evolution of the peaks for the selected genes over time'''

import argparse
import os
import sys
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;
import numpy as np
from collections import defaultdict



parser = argparse.ArgumentParser(description='Shows an evolution of the peaks for the selected genes over time');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated consensus regions, gff file");
parser.add_argument('--expression', nargs = '+', required = True, type = str, help = "Path to the gene expression, tsv format (differential/*.tsv) ");
parser.add_argument('--genes', nargs = '+', default=['ctaE', 'sigC', 'cydA'], required = False, type = str, help = "Path to the selected genes to show the binding evolution");
parser.add_argument('--table', nargs = '?', required = True, type = str, help = "Table which connects ncbi gene names and genesymbols, csv format");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();

table = {};


with open(args.table) as f:
    next(f);
    for l in f:
        a = l.strip().split(';')
        if(len(a)>2):
            table[a[2].split(' ')[0]] = a[1];
            
            
        
genes = [table[x] for x in args.genes];
selection = [0,2,3,4]
print("\t".join(['gene', '0h', '0.5h', '2h', '4h'])) 

chap = {}
for interval in BedTool(args.path):
    for symbol, gene in zip(args.genes, genes):
        if(gene in interval.attrs['start_gene']):
            a = np.array([1 if x == 'None' else float(x) for x in interval.attrs['maxcov'].split(",") ]);
            print("%s\t%s" % (symbol, "\t".join( [str(x) for x in a[selection] ])) ) ;
            chap[symbol] = a[selection]
            
print()
print(chap)

wt = defaultdict(list);
ko = defaultdict(list);

for efile in args.expression:
    with open(efile) as f:
        next(f);
        for l in f:
            for symbol, gene in zip(args.genes, genes):
                if(gene in l):
                    a = l.strip().split("\t")
                    wt[symbol].append(np.mean([float(x) for x in a[1].split(";")]));
                    ko[symbol].append(np.mean([float(x) for x in a[2].split(";")]));
print()
print(wt);
print(ko)

xdata1 = [0,1,2,3]
xdata2 = [0,1,3]

for symbol in args.genes:
    
    fig, ax1 = plt.subplots(figsize=(16,9))
    ax2 = ax1.twinx()

    ax1.set_xlabel('Time (h)', fontsize='x-large')
    ax1.set_ylabel('Max peak intensity', fontsize='x-large')
    ax2.set_ylabel('mRNA level (TPM)', fontsize='x-large')
    
    ax1.tick_params(axis='both', labelsize='x-large', top=False, right=False)
    ax2.tick_params(axis='both', labelsize='x-large', top=False, right=False)
    
    ax1.spines['top'].set_visible(False)
    #ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    #ax2.spines['right'].set_visible(False)
    plt.xticks(xdata1, ['0h', '0.5h', '2h', '4h'], fontsize= 'x-large')
    #ax1.set_xticklabels(['0h', '0.5h', '2h', '4h'], fontsize = 'x-large')

    ax1.plot(xdata1, chap[symbol], color = 'red', linestyle='dashed', linewidth=3, label='HrrA binding', marker='o', markersize = 10)
    ax2.plot(xdata2, wt[symbol], color = 'lightblue',linewidth=3, label='WT', marker='o', markersize = 10)
    ax2.plot(xdata2, ko[symbol], color = 'darkblue', linewidth=3, label='HrrA KO', marker='o', markersize = 10)
    fig.legend(loc=(0.15, 0.8), frameon=False, fontsize='x-large')




    #if(args.plot):
        #_format = args.plot.split(".")[-1]
        #plt.savefig(args.plot, format = _format)
    #else:
    plt.savefig(os.path.join(args.outdir, "total.%s.%s"  % (symbol, args.format)) , format = args.format)











































            
