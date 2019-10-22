#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Adjusts genomic coverage based on the provided control (dna expresssion) coverage'''

import argparse
import os
import sys
import numpy as np;
import pandas as pd;
from scipy.stats import pearsonr

from afbio.sequencetools import sliding_window
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Adjusts genomic coverage based on the provided control (dna expresssion) coverage');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the coverage file, bed3 format");
parser.add_argument('--control', nargs = '?', required=True, type = str, help = "Path to the control file, bed3 format");
parser.add_argument('--smooth', nargs = '?', default=0, type = int, help = "Sliding window half-length used to smooth the control genomic coverage, default: 0");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output plot directory")
args = parser.parse_args();

def smooth_coverage(coverage, flen):
    length = 2*flen+1
    res = [np.mean(x) for x in sliding_window(coverage, length)]
    res = [res[0]]*flen + res + [res[-1]]*flen
    return np.array(res);

####Read Coverage
data = pd.read_csv(args.path, sep="\t" , names = ["chr", "position", "coverage"])
coverage = data.coverage.values



control_coverage = pd.read_csv(args.control, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values
#control_coverage = smooth_coverage(control_coverage, args.smooth)
control_coverage += 1;
control_coverage = control_coverage/np.mean(control_coverage)

#data.coverage = coverage/control_coverage;
#for l in data.itertuples():
    #print("%s\t%s\t%1.1f" % (l.chr, l.position, l.coverage));



data.coverage = coverage/control_coverage;
for l in data.itertuples():
    print("%s\t%s\t%1.1f" % (l.chr, l.position, l.coverage));
    
 
 
### PLOTTING ### 
fontsize=20
linewidth = 5 

###Plot coverage vs genomic control correlation
size=500
correlations = []
for r, c in zip(np.array_split(coverage, size), np.array_split(control_coverage, size)):
    correlations.append(pearsonr(r, c)[0])
    
fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.95])
ax.set_xlabel("Pearson R", fontsize=fontsize)
ax.set_ylabel('Density', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=fontsize)
ax.tick_params(axis='both', which='minor', labelsize=fontsize)

ax.spines['top'].set_visible(False);
ax.spines['right'].set_visible(False)
ax.hist(correlations, density=True);
plt.savefig(os.path.join(args.outdir, "control_correlation.%s"  %  args.format) , format = args.format)
    
    
    
###Plot coverage vs convolution
size=200
smoothed_coverage = [];
for c in np.array_split(control_coverage, size):
    smoothed_coverage.extend([np.mean(c)]*len(c))




ax.plot(smoothed_coverage, 'b-', label='genomic control')
ax.set_xlabel("position (nt)", fontsize=fontsize)
ax.set_ylabel('genomic control', fontsize=fontsize)
ax2 = ax.twinx()
ax2.plot(coverage, 'r-', label="CHAP coverage")
ax2.set_ylabel('coverage', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=fontsize)
ax.tick_params(axis='both', which='minor', labelsize=fontsize)
ax2.tick_params(axis='both', which='major', labelsize=fontsize)
ax2.tick_params(axis='both', which='minor', labelsize=fontsize)

ax.spines['top'].set_visible(False) 
ax2.spines['top'].set_visible(False) 

fig.legend(frameon = False, fontsize=fontsize, loc=(0.65, 0.75))
plt.savefig(os.path.join(args.outdir, "control_coverage.%s"  %  args.format) , format = args.format)
    #plt.show()
