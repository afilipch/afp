#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Finds and explores differentially expressed genes'''

import argparse
import os
import sys
from collections import defaultdict
from itertools import combinations
from multiprocessing.dummy import Pool

import numpy as np;
import matplotlib.pyplot as plt;
from scipy.stats import fisher_exact, variation, pearsonr

from afbio.LRGFDR import lrg

#from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Finds and explores differentially expressed genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the gene expression file, tsv format (compile_expression.py output)");
parser.add_argument('--minexpr', nargs = '?', default=0.001, type = float, help = "Minimum required expression for the differential analysis");
parser.add_argument('--fdr', nargs = '?', default=0.05, type = float, help = "Threshold for the local empirical false discovery rate");
parser.add_argument('--plot', nargs = '?', required=True, type = str, help = "Output directory for the plots");
parser.add_argument('--labelnames', nargs = '?' , type = str, help = "Path to the file with gene names to be labeled on loglog plot");
args = parser.parse_args();


###########################################################################################################################################################
###PARSE INPUT FILES

basename  = os.path.basename(args.path).split(".")[0]

genenames = [];
expressions = [];
misfits = []

with open(args.path) as f:
    header = next(f).strip().split("\t")
    ratiostr = "%s/%s" % (header[2], header[1])
    for l in f:
        a = l.strip().split("\t")
        s1 = [float(x) for x in a[1].split(";")]
        s2 = [float(x) for x in a[2].split(";")]
        if(all([x>args.minexpr for x in s1+s2])):
            genenames.append(a[0]);
            expressions.append((s1, s2));
        else:
            misfits.append((a[0], s1, s2))

labelnames = [];
if(args.labelnames):
    with open(args.labelnames) as f:
        for l in f:
            labelnames.append(l.strip());
    


###########################################################################################################################################################
###Assign and explore intra- and inter-sample VARIATIONS

def get_variation(expression):
    return variation(expression[0]), variation(expression[1]), variation(expression[0] + expression[1])

variations = np.array([get_variation(x) for x in expressions])
variations = variations.transpose()
 
#Plot variations
plt.style.use('seaborn-deep')
bins = np.linspace(0,1,51)

plt.hist([variations[0], variations[2]], bins, label=['intra', 'inter'])
plt.legend(loc='upper right')
plt.savefig(os.path.join(args.plot, '%s.variations_1.png' % basename), format = 'png')
plt.clf()

plt.hist([variations[1], variations[2]], bins, label=['intra', 'inter'])
plt.legend(loc='upper right')
plt.savefig(os.path.join(args.plot, '%s.variations_2.png'  % basename), format = 'png')
plt.clf()

###########################################################################################################################################################
### Assign and explore intra- and inter-sample FOLD CHANGES

def localfold(a1, a2):
    return (a2-a1)/(a2+a1)

def get_fold(expression):
    return localfold(*expression[0]), localfold(*expression[1]), localfold(sum(expression[0]), sum(expression[1]))

def norm_fold_2_log_fold(fold):
    return np.log2( (fold+1)/(1-fold) )

folds = np.array([get_fold(x) for x in expressions])
folds = folds.transpose()
 
#Plot folds
plt.style.use('seaborn-deep')
bins = np.linspace(-1,1,51)

plt.hist([folds[0], folds[2]], bins, label=['intra', 'inter'])
plt.legend(loc='upper right')
plt.savefig(os.path.join(args.plot, '%s.folds_1.png' % basename), format = 'png')
plt.clf()

plt.hist([folds[1], folds[2]], bins, label=['intra', 'inter'])
plt.legend(loc='upper right')
plt.savefig(os.path.join(args.plot, '%s.folds_2.png' % basename), format = 'png')
plt.clf()


###########################################################################################################################################################
#GET COMPOSITE THRESHOLDS based on fold change and expression


noise_fold = [max( abs(x[0]), abs(x[1]) ) for x in zip(folds[0], folds[1])]
signal_fold = [abs(x) for x in folds[2]]
logfolds = [norm_fold_2_log_fold(x) for x in folds[2]]

common_expr = [np.log2( (max(x[0]) + max(x[1]))/2.0 + 1)  for x in expressions]
noise_total =  list(zip(noise_fold, common_expr))
signal_total = list(zip(signal_fold, common_expr))

def get_fold_cutoff(signal, noise, fdr):
    threshold = 1.01;
    for t in np.linspace(0,1,1001):
        fs = len([x for x in signal if x > t])
        fn = len([x for x in noise if x > t])
        if(not fs):
            return threshold, fs, fn
        p = fn/(fs+fn);
        if(p<fdr):
            threshold = t;
            return threshold, fs, fn
    return threshold, fs, fn;


levels = [np.percentile(common_expr, x) for x in np.linspace(0,100,21)]
levels[0] -= 0.001
thresholds = [];
for l1, l2 in zip(levels, levels[1:]):
    signal = [x[0] for x in signal_total if x[1]>l1 and x[1]<=l2]
    noise = [x[0] for x in noise_total if x[1]>l1 and x[1]<=l2]
    threshold, ps, pn = get_fold_cutoff(signal, noise, args.fdr)
    sys.stderr.write("%.1f\t%.1f\t%.2f\t%d\t%d\t%d\t%d\n" % (l1, l2, threshold, len(signal), ps, len(noise), pn) )
    thresholds.append((l1, l2, threshold))
    
lboundary = levels[0];
for l1, l2, t in thresholds:
    if(t<1):
        lboundary = l1;
        break;
    
#print(lboundary)
#sys.exit()



#upper_signal = [x[0] for x in signal_total if x[1]>lboundary]
#upper_noise = [x[0] for x in noise_total if x[1]>lboundary]
#basal_level = get_fold_cutoff(upper_signal, upper_noise, args.fdr)[0]

def trline(expr, lboundary, factor, basal_level):
    if(expr>lboundary):
        return 2**(-factor*expr)+basal_level
    else:
        return 1;

def estimate_fdr_for_trline(signal, noise, factor, basal_level):
    fs = len([x for x in signal if x[0] > trline(x[1], lboundary, factor, basal_level)])
    fn = len([x for x in noise if x[0] > trline(x[1], lboundary, factor, basal_level)  ])
    if(fn):
        return fs, fn, fn/(fs+fn)
    else:
        return fs, 0, 0;
    
    

basal_level = list(sorted([x[2] for x in thresholds]))[0]
roc = [];
for factor in np.linspace(0, 2, 41):
    fs, fn, fdr= estimate_fdr_for_trline(signal_total, noise_total, factor, basal_level)
    #print("%.2f\t%d\t%d\t%.3f" % (factor, fs, fn, fdr))
    roc.append((fs/len(signal_total), fdr, factor))

#rough strategy so far
sensitivity, empirical_fdr, bestfactor = max([x for x in roc if x[1]<args.fdr], key=lambda x: x[0])
sys.stderr.write("\nsensitivity: %.3f\nFDR: %.3f\nselected factor: %.3f\nselected left boundary: %.3f\nselected basal level: %.3f" % (sensitivity, empirical_fdr, bestfactor, lboundary, basal_level));
#sys.exit();



### step-wise approach
#levels[-1] += 0.01
#thresholds = [];
#for l1, l2 in zip(levels, levels[1:]):
    #signal = [x[0] for x in signal_total if x[1]>=l1 and x[1]<l2]
    #noise = [x[0] for x in noise_total if x[1]>=l1 and x[1]<l2]
    #threshold, ps, pn = get_fold_cutoff(signal, noise, args.fdr)
    #sys.stderr.write("%.1f\t%.1f\t%.2f\t%d\t%d\t%d\t%d\n" % (l1, l2, threshold, len(signal), ps, len(noise), pn) )
    #thresholds.append((l1, l2, threshold))
    




###########################################################################################################################################################
### Output the results

def check_diff(normed_fold, common):
    tr = trline(common, lboundary, bestfactor, basal_level) 
    return int(normed_fold>tr)

scatter = []
textlabels = [];

print("\t".join(header + ["log2(%s)" % ratiostr, "normed_fold", "diff"]))
for name, expression, logfold, normed_fold, common_expr in zip(genenames, expressions, logfolds, folds[2], common_expr):
    e1 = ";".join([str(x) for x in expression[0]])
    e2 = ";".join([str(x) for x in expression[1]])
    diff = check_diff(abs(normed_fold), common_expr)
    print("%s\t%s\t%s\t%1.3f\t%1.3f\t%d" % (name, e1, e2, logfold, normed_fold, diff))
    scatter.append((common_expr, logfold, diff))
    
for name, wt_values, ko_values in misfits:
    e1 = ";".join([str(x) for x in wt_values])
    e2 = ";".join([str(x) for x in ko_values])
    wt = np.mean(wt_values)
    ko = np.mean(ko_values)
    if(wt+ko):
        localfold = (ko-wt)/(wt+ko)
    else:
        localfold = 0;
    if(wt and ko):
        logfold = np.log2(ko/wt)
    elif(ko):
        logfold = 999999
    elif(wt):
        logfold = -999999
    else:
        logfold = 0;
    print("%s\t%s\t%s\t%1.3f\t%1.3f\t%d" % (name, e1, e2, logfold, localfold, 0))
    
        
    

    if(name in labelnames):
        textlabels.append((name, common_expr, logfold))







    
###generate loglog plot
scatter1 = np.array([(x[0], x[1]) for x in scatter if not x[-1]]);
scatter2 = np.array([(x[0], x[1]) for x in scatter if x[-1]]);
sys.stderr.write("\n%d genes are differentially expressed\n" % scatter2.shape[0])

fig, ax = plt.subplots(figsize=(16, 9))
plt.plot(scatter1[:,0], scatter1[:,1], marker = 'o', markerfacecolor='lightblue', markeredgecolor='lightblue', linestyle = "None", markersize = 4, label = 'non-differential');
plt.plot(scatter2[:,0], scatter2[:,1], marker = 'o', markerfacecolor='coral', markeredgecolor='coral', linestyle = "None", markersize = 4, label = "differential");
plt.ylabel("log2(%s)" % ratiostr);
plt.xlabel("log2(TPM)")
plt.legend(frameon=False, fontsize='xx-large')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.title("Expression vs Fold Change")

for text, x, y in textlabels:
    ax.text(x, y, text);

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize('xx-large')


plt.savefig(os.path.join(args.plot, '%s.scatter.eps' % basename), format = 'eps')
plt.clf()


###Generate scatter expression plot
yvals = [np.log2(sum(x[0])+1) for x in expressions]
xvals = [np.log2(sum(x[1])+1) for x in expressions]
fig, ax = plt.subplots(figsize=(16, 9))
plt.plot(xvals, yvals, marker = 'o', markerfacecolor='lightblue', markeredgecolor='lightblue', linestyle = "None", markersize = 4);
plt.xlabel("log2(%s)" % header[1]);
plt.ylabel("log2(%s)" % header[2]);
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.title("%s vs %s" % tuple(header[1:3]))
pcorr_log = pearsonr(xvals, yvals)[0]
pcorr = pearsonr([sum(x[0]) for x in expressions], [sum(x[1]) for x in expressions])[0]
sys.stderr.write("\nLog2 correlation coefficent between %s and %s is equal %.3f\n" % (header[1], header[2], pcorr_log))
sys.stderr.write("\nCorrelation coefficent between %s and %s is equal %.3f\n" % (header[1], header[2], pcorr))
plt.text(0.1, 0.8, 'r(log)=%.2f\nr=%.2f' % (pcorr_log, pcorr), transform=ax.transAxes, fontsize = 'xx-large')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize('xx-large')
plt.savefig(os.path.join(args.plot, '%s.loglog.eps' % basename), format = 'eps')
plt.clf()




              













