#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explore evolution of the peaks over time'''

import argparse
import os
import sys
from itertools import combinations;



import numpy as np;
from scipy.stats import pearsonr, normaltest, variation
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;



parser = argparse.ArgumentParser(description='Explore evolution of the peaks over time');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the consensus regions, gff file");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
args = parser.parse_args();

intensities_list = [];
for interval in BedTool(args.path):
    intensity = [float(x) for x in interval.attrs['topcoverage'].split(",")]
    if(all(intensity)):
        intensities_list.append(intensity)
 
intensities_list = np.array(intensities_list)
low = np.percentile(intensities_list, 5, axis=1)
high = np.percentile(intensities_list, 95, axis=1)
width = high - low
means = np.mean(intensities_list, axis=1)

all_pairs = list(zip(width, means))
normed_width = [x[0]/x[1] for x in all_pairs]
#print(variation(normed_width), variation(width))
print(pearsonr(width, means))
print(pearsonr(normed_width, means))
print(pearsonr(width, np.log(means)))



strong_pairs = [x for x in all_pairs if x[1]>=10]
strong_width = [x[0] for x in strong_pairs]
strong_normed_width = [x[0]/x[1] for x in strong_pairs]
#print(pearsonr([x[0] for x in strong_pairs], [x[1] for x in strong_pairs]))
weak_pairs = [x for x in all_pairs if x[1]<10]
weak_normed_width = [x[0]/x[1] for x in weak_pairs]


#normality = normaltest(intensities_list, axis=1).pvalue
#normality.sort()
#print(["%1.5f" % x for x in normality])
    

print(len(strong_width))

### DRAWING SECTION ###

fontsize = 20
fig, ax = plt.subplots(figsize=(16,9))

ax.set_xlabel('peak intensity', fontsize=fontsize)
ax.set_ylabel("Width of 0.95 confidence interval", fontsize=fontsize)    
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.plot(means, width, 'b.')   

plt.savefig(os.path.join(args.outdir, "confidence_intensity_correlation.%s"  %  args.format ) , format = args.format)
plt.clf()
plt.close




fig, ax = plt.subplots(figsize=(16,9))
ax.set_ylabel("Normalized width of 0.95 confidence interval", fontsize=fontsize)    
ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.boxplot([strong_normed_width, weak_normed_width], notch = True, showfliers=True, labels = ['strong', 'weak'])   
#plt.xticks([], [])
plt.savefig(os.path.join(args.outdir, "boxplot.%s"  %  args.format ) , format = args.format)
plt.clf()
plt.close




#fig, ax = plt.subplots(figsize=(16,9))
#ax.set_ylabel("Normalized width of 0.95 confidence interval", fontsize=fontsize)    
#ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.boxplot(weak_normed_width, notch = True, showfliers=True)   
#plt.xticks([], [])
#plt.savefig(os.path.join(args.outdir, "weak_boxplot.%s"  %  args.format ) , format = args.format)
#plt.clf()
#plt.close
