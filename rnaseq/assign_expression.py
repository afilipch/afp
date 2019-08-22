#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Assigns expression (TPM, RPKM/FPKM) to the genes based on provided mappings (sorted bed file)'''

import argparse
import os
import sys
#from collections import defaultdict
import math

import numpy as np;
from scipy.stats import variation 
import matplotlib.pyplot as plt;
from scipy.stats.stats import pearsonr;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Assigns expression (TPM, RPKM/FPKM) to the genes based on provided mappings (sorted bed file)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the mapped reads, sorted bed format");
parser.add_argument('--genes', nargs = '?', required=True, type = str, help = "Path to the genes, sorted gff format");
parser.add_argument('--plot', nargs = '?', type = str, help = "Output destination for the statistics plots");
args = parser.parse_args();

def assign_coverage(intervals):
    coverage = np.array([float(x.fields[10]) for x in intervals])
    if(not any(coverage)):
        return[0, 0, len(intervals[0])/1000, 0, 'None']
    rpk = coverage.mean()*1000;
    return [rpk, sum(coverage), len(coverage)/1000, coverage.mean(), variation(coverage)]


mapping = BedTool(args.path)
annotation = BedTool(args.genes)
name2stat = {};
curname = '';
curints = [];
for cov in annotation.coverage(b=mapping, F=0.2, s=True, sorted=True, d=True):
    if(cov.name == curname):
        curints.append(cov);
    else:
        if(curints):
            name2stat[curname] = assign_coverage(curints);
        curints = [cov];
        curname = cov.name
else:
    name2stat[curname] = assign_coverage(curints);

    
### Normalize to RPKM
normfactor = sum([x[1] for x in name2stat.values() if x[1]])/1000000
for l in name2stat.values():
    l[1] = (l[1]/normfactor)/l[2]

    
### Normalize to TPM
normfactor = sum([x[0] for x in name2stat.values() if x[0]])/1000000
for l in name2stat.values():
    l[0] = l[0]/normfactor


for gene in annotation:
    stat = name2stat.get(gene.name, [0, 0, len(gene)/1000, 0, 'None']);
    gene.attrs['tpm'] = "%.2f" % stat[0]
    gene.attrs['rpkm'] = "%.2f" % stat[1]
    gene.attrs['covmean'] = "%.2f" % stat[3]
    sys.stdout.write(str(gene))
    
    
sys.stderr.write("\nTotal number of genes: %d\nNumber of expressed genes: %d\n\n" % (len(annotation), len([x for x in name2stat.values() if x[3]>=1])) )
    
    
if(args.plot):
    plotpath = os.path.join(args.plot, ".".join(os.path.basename(args.path).split(".")[:-1]))
    exppath = plotpath  + ".tpm_vs_rpkm.png"
    varpath = plotpath  + ".variation.png"
    tpmpath = plotpath  + ".tpm_distribution.png"
    
    plt.figure(figsize=(16,12))
    ax = plt.axes()
    statistics = name2stat.values();
    tpms = [x[0] for x in statistics if x[0]]
    rpkms = [x[1] for x in statistics if x[0]]
    r = pearsonr(tpms, rpkms)[0]
    plt.plot(rpkms, tpms, 'ro', color = 'darkblue');
    plt.xlabel('tpm')
    plt.ylabel('rpkm')
    plt.text(0.1, 0.8, "Pearson: %.2f%%" % (r*100), transform=ax.transAxes);
    plt.savefig(exppath, format='png');
    plt.clf();
    
    variations = [x[4] for x in statistics if x[0]]
    #sys.stderr.write("%s\n" % variations)
    plt.hist(variations, color='coral', bins=50, density=True)
    plt.xlabel('coverage coefficient of variation')
    plt.ylabel('fraction of genes')
    plt.savefig(varpath, format='png');
    plt.clf();
    
    #plt.xscale('log', basex=10)
    #ax = plt.axes()
    #ax = plt.axes()
    stub1, bins, stub2 = plt.hist([math.log(x,10) for x in tpms], color='coral', bins=50, density=True)
    plt.xlabel('log10 TPM')
    plt.ylabel('fraction of genes')
    #sys.stderr.write("%s\n" % bins)
    xticklabels = [int(10**x) for x in np.linspace(0, int(max(bins)), int(max(bins))+1) ]
    #sys.stderr.write("%s\n" % xticklabels)
    plt.xticks([math.log(x, 10) for x in xticklabels], [str(x) for x in xticklabels])
    plt.savefig(tpmpath, format='png');
    
    
    
    









