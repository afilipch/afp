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



parser = argparse.ArgumentParser(description='Detects peaks from the covergage track using a kernel-based convolution');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the coverage file");
parser.add_argument('--kernel', nargs = '?', choices = ['ndg', 'square'], default='ndg', type = str, help = "Kernel to use for the fitering");
parser.add_argument('--peakwidth', nargs = '?', default=0, type = int, help = "Expected width of the peaks, if not set peakwidth is estimated automatically");
parser.add_argument('--widthfactor', nargs = '?', default=1, type = float, help = "Multiplier for the eastimated peakwidth, so peakwidth=estimated_peakwidth*widthfactor");
parser.add_argument('--meanmult', nargs = '?', default=6.0, type = float, help = "Minimum required peak height (denoted in median coverage) for the peakwidth estimation");
parser.add_argument('--threads', nargs = '?', default=8, type = int, help = "Number of threads");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");
parser.add_argument('--convolution', nargs = '?', default='', type = str, help = "If set, the convolution track is written to the provided file");
args = parser.parse_args();

###Read Coverage
coverage_dict = coverage2dict(args.path);
#print( [x[:100] for x in coverage_dict.values()] )
coverage = []
for v in coverage_dict.values():
    coverage.extend(v)

###Output basic coverage statistics
sys.stderr.write("###detect_peaks\n")
sys.stderr.write("Median coverage:\t%1.2f\nMean coverage:\t%1.2f\nCoverage standard deviation:\t%1.2f\nMax coverage:\t%1.2f\n\n" % (np.median(coverage), np.mean(coverage), np.std(coverage), max(coverage)))



###Determine bandwidth of the potential peaks
#valid = True
if(args.peakwidth):
    bandwidth = args.peakwidth
    sys.stderr.write("Bandwidth was manually set up to:\t%1.2f\n\n" % bandwidth)
else:
    rawbandwidth, numpeaks = estimate_bandwidth(coverage, args.meanmult)
    if(numpeaks < 10):
        ##sys.stdout.write('');
        bandwidth = 120
        sys.stderr.write("\nWARNING! Only %d peak(s) detected. Bandwidth is set to the default value: %d\n\n" % (numpeaks, bandwidth))
        sys.stderr.write("Raw bandwith:\t%1.2f\nAdjusted bandwith:\t%1.2f\nRaw peaks detected:\t%1.2f\n\n" % (bandwidth, bandwidth, numpeaks))
        #valid = False
    else:
        bandwidth = rawbandwidth*args.widthfactor
        sys.stderr.write("Raw bandwith:\t%1.2f\nAdjusted bandwith:\t%1.2f\nRaw peaks detected:\t%1.2f\n\n" % (rawbandwidth, bandwidth, numpeaks))

#if(valid):

###Convolute coverage with a given kernel
exec("from afbio.filters import %s as kernelfunc" % args.kernel)
kernel = kernelfunc(truncate = 4.0);


convolution_list = [];
total_peaks= 0;
for chrom, ch_cov in coverage_dict.items():
    normed_cov = ch_cov/np.mean(ch_cov)
    convolution = np.array(convolute(normed_cov, kernel, bandwidth, threads=args.threads))
    convolution_list.append((chrom, convolution));
    peaks = detect_peaks(convolution);
    total_peaks += len(peaks)
    for pk in peaks:
        sys.stdout.write("%s\t%d\t%d\t%d\t%1.3f\t%s\n" % (chrom, pk[0], pk[2], pk[1], pk[3], '+'));


###Output basic statistics of the detected peaks
sys.stderr.write("\nRaw peaks detected with a kernel:\t%d\n\n" % len(peaks))

###Output convolution track
if(args.convolution):
    with open(args.convolution, 'w') as f:
        for chrom, convolution in convolution_list:
            for c, el in enumerate(convolution, start=1):
                f.write("%s\t%d\t%d\n" % (chrom, c, el)); 

###Plot coverage vs convolution
if(args.plot):
    import matplotlib.pyplot as plt
    fontsize = 24
    fig, ax1 = plt.subplots(figsize=(16,9))
    plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.95])
    

    ax1.plot(coverage, 'b-')
    ax1.set_xlabel("position (nt)", fontsize=fontsize)
    ax1.set_ylabel('coverage', color='b', fontsize=fontsize)
    ax1.tick_params('y', colors='b')
    
    convolution = [max(0,x) for x in convolution]
    ax2 = ax1.twinx()
    ax2.plot(convolution, 'r-')
    ax2.set_ylabel("convolution", color='r', fontsize=fontsize)
    ax2.tick_params('y', colors='r')
    
    #ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False) 
    ax2.spines['top'].set_visible(False) 
    ax1.tick_params(axis='both', which='major', labelsize=fontsize)
    ax1.tick_params(axis='both', which='minor', labelsize=fontsize)
    ax2.tick_params(axis='both', which='major', labelsize=fontsize)
    ax2.tick_params(axis='both', which='minor', labelsize=fontsize)


    
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format)
    #plt.show()

    
        
