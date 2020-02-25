#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Detects peaks from the covergage track using a kernel-based convolution'''

import argparse
import os
import sys
#import scipy
import numpy as np;
import pandas as pd;
from afbio.peaks import estimate_bandwidth, convolute, detect_peaks
from afbio.sequencetools import coverage2dict
from multiprocessing.dummy import Pool

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch



parser = argparse.ArgumentParser(description='Detects peaks from the covergage track using a kernel-based convolution');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the coverage file");
parser.add_argument('--kernel', nargs = '?', choices = ['ndg', 'square'], default='ndg', type = str, help = "Kernel to use for the fitering");
parser.add_argument('--outdir', nargs = '?', type = str, help = "Path to the output plot folder");
parser.add_argument('--format', default = 'png', nargs = '?', type = str, help = "Plot format");
args = parser.parse_args();


###Read Coverage
coverage_dict = coverage2dict(args.path);
coverage = coverage_dict['NC_003450.3'][1661900:1667000]
bandwidth = 120



#if(valid):

###Convolute coverage with a given kernel
exec("from afbio.filters import %s as kernelfunc" % args.kernel)
kernel = kernelfunc(truncate = 4.0);

convolution = np.array(convolute(coverage, kernel, bandwidth, threads=2))
peaks = detect_peaks(convolution);



#########################################################################################################################
### Drawing section



colors = ["deepskyblue", "darkmagenta", "limegreen"]
fontsize = 24

### Draw convolution ###
fig, ax1 = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.95])


ax1.plot(coverage, color = colors[0])
ax1.set_xlabel("position (nt)", fontsize=fontsize)
ax1.set_ylabel('coverage', color=colors[0], fontsize=fontsize)
ax1.tick_params('y', colors=colors[0])

convolution = np.array([max(0,x) for x in convolution])
convolution[0:400] *= 0.1
convolution[-400:] *= 0.1;
ax2 = ax1.twinx()
ax2.plot(convolution, color = colors[1])
ax2.set_ylabel("convolution", color=colors[1], fontsize=fontsize)
ax2.tick_params('y', colors=colors[1])

#ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False) 
ax2.spines['top'].set_visible(False) 
ax1.tick_params(axis='both', labelsize=fontsize)
ax2.tick_params(axis='both', labelsize=fontsize)

plt.savefig(os.path.join(args.outdir, "convolution.%s" %  args.format), format = args.format);
plt.clf()
plt.close()


### Draw kernel ###

fig, ax = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.95])

#print(len(kernel))
norma = 0.5*max(coverage)/max(kernel)
y_kernel = np.array([kernel[int(x)] for x in np.linspace(0, len(kernel), 1001)[:-1]])*norma + norma


ax.plot(coverage, color = colors[0])
ax.plot(y_kernel, color = colors[1])
ax.arrow(len(y_kernel), (max(y_kernel) + min(y_kernel))/2, len(y_kernel)/2, 0, lw = 6, color = colors[1], head_width=12, head_length=50)
ax.text(0.5, 1, "Convolution with second-order gaussian kernel enhances peaks and supresses noise", horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=fontsize*0.8)


ax.set_xlabel("position (nt)", fontsize=fontsize)
ax.set_ylabel('coverage', fontsize=fontsize)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)  
ax.tick_params(axis='both', labelsize=fontsize)

plt.savefig(os.path.join(args.outdir, "kernel.%s" %  args.format), format = args.format);
plt.clf()
plt.close()


### Draw filtering
THRESHOLD = 1200;

fig, ax1 = plt.subplots(figsize=(16,9))
plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.95])

ax1.plot(coverage, color = colors[0])
ax1.set_xlabel("position (nt)", fontsize=fontsize)
ax1.set_ylabel('coverage', color=colors[0], fontsize=fontsize)
ax1.tick_params('y', colors=colors[0])

convolution = np.array([max(0,x) for x in convolution])
convolution[0:400] *= 0.1
convolution[-400:] *= 0.1;
ax2 = ax1.twinx()

start = 0;
for peak in peaks:
    if(convolution[peak[1]] > THRESHOLD):
        ax2.plot(range(start, peak[0]), convolution[start:peak[0]], color = colors[1])
        ax2.plot(range(peak[0], peak[2]), convolution[peak[0]:peak[2]], color = colors[2])
        start = peak[2];
else:
    ax2.plot(range(start, len(convolution)), convolution[start:], color = colors[1])

ax2.axhline(y=THRESHOLD, color = colors[2], lw = 4, ls = '--')
ax2.text(len(convolution)*0.53, THRESHOLD*1.1, "threshold with a controlled false discovery rate", fontsize=fontsize*0.7)

ax2.set_ylabel("convolution", color=colors[1], fontsize=fontsize)
ax2.tick_params('y', colors=colors[1])

#ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False) 
ax2.spines['top'].set_visible(False) 
ax1.tick_params(axis='both', labelsize=fontsize)
ax2.tick_params(axis='both', labelsize=fontsize)



        
    #ax2.arrow(peak[1], top-40, 0, 20, lw = 6, color = 'purple', head_width=12, head_length=20)

plt.savefig(os.path.join(args.outdir, "filtering.%s" %  args.format), format = args.format);
plt.clf()
plt.close()





    
        
