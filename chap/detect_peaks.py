#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Detects peaks from the covergage track using a kernel-based convolution'''

import argparse
import os
import sys
#import scipy
import numpy as np;
import pandas as pd;
from afbio.peaks import estimate_bandwidth, convolute, detect_peaks
from multiprocessing.dummy import Pool


parser = argparse.ArgumentParser(description='Detects peaks from the covergage track using a kernel-based convolution');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the coverage file");
parser.add_argument('--kernel', nargs = '?', choices = ['ndg', 'square'], default='ndg', type = str, help = "Kernel to use for the fitering");
parser.add_argument('--peakwidth', nargs = '?', default=0, type = int, help = "Expected width of the peaks, if not set peakwidth is estimated automatically");
parser.add_argument('--widthfactor', nargs = '?', default=1, type = float, help = "Multiplier for the eastimated peakwidth, so peakwidth=estimated_peakwidth*widthfactor");
parser.add_argument('--meanmult', nargs = '?', default=6.0, type = float, help = "Minimum required peak height (denoted in median coverage) for the peakwidth estimation");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");
parser.add_argument('--convolution', nargs = '?', default='', type = str, help = "If set, the convolution track is written to the provided file");
args = parser.parse_args();

###Read Coverage
coverage = pd.read_csv(args.path[0], sep="\t" , names = ["chr", "position", "coverage"]).coverage.values
for path in args.path[1:]:
    coverage += pd.read_csv(path, sep="\t" , names = ["chr", "position", "coverage"]).coverage.values

###Output basic coverage statistics
sys.stderr.write("Median coverage:\t%1.2f\nMean coverage:\t%1.2f\nCoverage standard deviation:\t%1.2f\nMax coverage:\t%1.2f\n\n" % (np.median(coverage), np.mean(coverage), np.std(coverage), max(coverage)))

###Determine bandwidth of the potential peaks
if(args.peakwidth):
    bandwidth = args.peakwidth
    sys.stderr.write("Bandwidth was manually set up to:\t%1.2f\n\n" % bandwidth)
else:
    rawbandwidth, numpeaks = estimate_bandwidth(coverage, args.meanmult)
    bandwidth = rawbandwidth*args.widthfactor
    sys.stderr.write("Raw bandwith:\t%1.2f\nAdjusted bandwith:\t%1.2f\nRaw peaks detected:\t%1.2f\n\n" % (rawbandwidth, bandwidth, numpeaks))


###Convolute coverage with a given kernel
exec("from afbio.filters import %s as kernelfunc" % args.kernel)
kernel = kernelfunc(truncate = 4.0);
#sys.stderr.write("bb\n")
convolution = np.array(convolute(coverage, kernel, bandwidth, threads=8))
#sys.stderr.write("%s\n" % convolution)


###Detect peaks for convolved coverage
peaks = detect_peaks(convolution);
    
###Outut the detected peaks
#sys.stdout.write("start\ttop\tend\tscore\n");
for pk in peaks:
    sys.stdout.write("%s\t%d\t%d\t%d\t%1.2f\t%s\n" % ('chr1', pk[0], pk[2], pk[1], pk[3], '+'));
 
 
 
###Output basic statistics of the detected peaks
sys.stderr.write("\nRaw peaks detected with a kernel:\t%d\n\n" % len(peaks))
 
###Output convolution track
if(args.convolution):
    with open(args.convolution, 'w') as f:
        for c, el in enumerate(convolution):
            f.write("chr1\t%d\t%d\n" % (c+1, el)); 
 
###Plot coverage vs convolution
if(args.plot):
    import matplotlib.pyplot as plt

    fig, ax1 = plt.subplots()

    ax1.plot(coverage, 'b-')
    ax1.set_xlabel("postion (nt)")
    ax1.set_ylabel('coverage', color='b')
    ax1.tick_params('y', colors='b')
    
    convolution = [max(0,x) for x in convolution]
    ax2 = ax1.twinx()
    ax2.plot(convolution, 'r-')
    ax2.set_ylabel("convolution", color='r')
    ax2.tick_params('y', colors='r')
    
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False) 

    fig.tight_layout()
    
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format)
    
    
        