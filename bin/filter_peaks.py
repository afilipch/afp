#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Filters the detected peaks based on the distribution of their scores'''

import argparse
import os
import sys
import scipy
import numpy as np;
import pandas as pd;
from scipy.optimize import curve_fit;
import matplotlib.pyplot as plt;


from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Filters the detected peaks based on the distribution of their scores');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the detected peaks");
parser.add_argument('--zscore', nargs = '?', default=3.0, type = float, help = "Z-score of the filtered peaks regarding to the fitted gaussian distribution");
parser.add_argument('--pcut', nargs = '?', default=99.0, type = float, help = "Percentile cuttof apllied to the peak scores prior fitting them to gaussian distribution");
parser.add_argument('--coverage', nargs = '?', default='', type = str, help = "Path to the coverage file, if provided peaks are filtered according to --minmedian option");
parser.add_argument('--minmedian', nargs = '?', default=3.0, type = float, help = "Minimal required height of the detected peaks, this filter is applied when --coverage is provided");

parser.add_argument('-ap', '--assignedpeaks', nargs = '?', default='', type = str, help = "If set, all the peaks with the assigned z-score are written to the given destination");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output plot (convolution scores distribution and filtering)");
args = parser.parse_args();


### Read the detected peaks
peaks = BedTool(args.path)
scores = [float(x.score) for x in peaks] 

### Output basic statistics of the read peaks
sys.stderr.write("%s\nfile is processed:\t%s\n\n" % ("_"*140, os.path.abspath(args.path)))
sys.stderr.write("median:\t%1.2f\nmean:\t%1.2f\nstd:\t%1.2f\n\n" % (np.median(scores), np.mean(scores), np.std(scores)))



### Get gaussian fit for the peak scores
pcut = np.percentile(scores, args.pcut)
cutscores = [x for x in scores if x < pcut]

suniq = len(set(scores))
nbins = min((1000, suniq))
hist, bin_edges = np.histogram(cutscores, density=True, bins = nbins)
bin_centers = (bin_edges[:-1] + bin_edges[1:])/2

## Define gauss function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

## p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
p0 = [max(hist), np.median(cutscores), np.median(cutscores)]

## Get fitting parameters
coeff, stub = curve_fit(gauss, bin_centers, hist, p0=p0, maxfev = 2000)
stub, fitmean, fitdev = coeff

## Print fitted parameters
sys.stderr.write('fitted mean:\t%1.2f\nfitted deviation:\t%1.2f\n\n' % (fitmean, fitdev))



### Filter the scores based on gaussian fitting
threshold = fitmean + args.zscore*fitdev
filtered = [x for x in peaks if float(x.score) > threshold]
scores_filtered = [float(x.score) for x in peaks if float(x.score) > threshold]

### Print basic statistics of the filtered peaks
sys.stderr.write("\nThreshold:\t%1.2f\nnum of raw peaks:\t%d\nnum of gauss-filtered peaks:\t%d\n\n" % (threshold, len(scores), len(filtered)))
sys.stderr.write("Filtered median:\t%1.2f\ninitial median:\t%1.2f\n\n" % (np.median(scores_filtered), np.median(scores)))

###Filter peaks based on the coverage properties
if(args.coverage):
    coverage = pd.read_csv(args.coverage, sep="\t" , names = ["chr", "postion", "coverage"]).coverage.values
    threshold = np.mean(coverage)*args.minmedian;
    before = len(filtered)
    filtered = [ x for x in filtered if max(coverage[x.start:x.end]) > threshold]
    sys.stderr.write("\nThreshold for min peak height:\t%1.2f\nnum of gauss-filtered peaks:\t%d\nnum of coverage-filtered peaks:\t%d\n\n" % (threshold, before, len(filtered)))  


 
### Print the filtered peaks
for peak in filtered:
    zscore = (float(peak.score) - fitmean)/fitdev;
    peak.score = "%1.1f" % zscore;
    sys.stdout.write(str(peak));
    

### Print z-values for all the peaks
if(args.assignedpeaks):
    with open(args.assignedpeaks, 'w') as f:
        for peak in peaks:
            peak.score = "%1.1f" % ((float(peak.score) - fitmean)/fitdev)
            f.write(str(peak))

### Plot the fitted gaussian curve vs scores distribution
if(args.plot):
    hist_fit = gauss(bin_centers, *coeff)
    plt.plot(bin_centers, hist, label='Test data')
    plt.plot(bin_centers, hist_fit, label='Fitted data')
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format)


