#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Correlates differential expression of the genes with the peak intensity on their promoters'''

import argparse
import os
import sys
#import copy
from collections import defaultdict, Counter
from itertools import product;

import numpy as np;
from scipy.stats import pearsonr
#import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

#from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Correlates differential expression of the genes with the peak intensity on their promoters');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated consensus regions of the time-series binding peaks");
parser.add_argument('--diff', nargs = 3, required=True, type = str, help = "Path to the file with genes\' differential expression");
parser.add_argument('--distance', nargs = '?', default=400, type = int, help = "Maximum allowed distance for the peak to the closest peak");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path to the plot directory");
args = parser.parse_args();

###READ differential gene expression
ldiff = len(args.diff)
def readdiff(path):
    tempd = {};
    with open(path) as f:
        next(f);
        for l in f:
            a = l.strip().split("\t");
            tempd[a[0]] = float(a[3])
    return tempd
    

names2diff = defaultdict(list);
for path in args.diff:
    for k, v in readdiff(path).items():
        names2diff[k].append(v);
names2diff = dict([x for x in names2diff.items() if len(x[1]) == ldiff]); 
        
print(names2diff)


###READ peak intensities        
        
names2intensities = {};        
for interval in BedTool(args.path):
    distance = int(interval.attrs['start_gene_distance'])
    if(distance < args.distance):
        genename = interval.attrs['start_gene']
        intensity = [float(x) if x !='None' else 0 for x in interval.attrs['maxcov'].split(",")]
        lint = len(intensity);
        names2intensities[genename] = intensity;


###CORRELATE peak intensities to differential gene expression
diff2timepoints = {0: '0h', 1: '0.5h', 2: '4h'}
int2timepoints = {0: 'pre', 1: '0h', 2: '0.5h', 3: '2h', 4: '4h', 5: '9h', 6: '24h'}

print("chap timepoint\trnaseq timepoint\tcorrelation for upregulated\tcorrelation for downregulated\tnum up\tnum down");


intensity2diff = {};
for i1, i2 in product(range(ldiff), range(lint)):
    up_int, down_int, up_diff, down_diff = [], [], [], []
    for genename, intensity in names2intensities.items():
        s_int = intensity[i2];
        if(s_int):
            diff = names2diff.get(genename, None);
            if(diff):
                s_diff = diff[i1];
                if(s_diff>0):
                    up_int.append(s_int);
                    up_diff.append(s_diff);
                elif(s_diff<0):
                    down_int.append(s_int);
                    down_diff.append(abs(s_diff));
    
    if(len(up_int)>5):
        up_rcoeff = pearsonr(up_int, up_diff)[0];
    else:
        up_rcoeff = 0;
    if(len(down_int)>5):
        down_rcoeff = pearsonr(down_int, down_diff)[0];
    else:
        down_rcoeff = 0;
        
    intensity2diff[(i2, i1)] = up_rcoeff, down_rcoeff
    print("%s\t%s\t%.3f\t%.3f\t%d\t%d" % (int2timepoints[i2], diff2timepoints[i1], up_rcoeff, down_rcoeff, len(up_int), len(down_diff)));
    
    
    
###PLOT 
selection = [(0,0), (0,1), (2,1), (0,2), (2,2), (4,2)]
upvals = [intensity2diff[x][0] for x in selection];
downvals = [intensity2diff[x][1] for x in selection];

ind = np.arange(len(selection))  # the x locations for the groups
width = 0.35       # the width of the bars

fig, ax = plt.subplots(figsize = (16, 9))
plt.tight_layout(rect=[0.04, 0.1, 1, 1])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
rects1 = ax.bar(ind, upvals, width, color='darkblue')
rects2 = ax.bar(ind + width, downvals, width, color='lightblue')

# add some text for labels, title and axes ticks
ax.set_ylabel('R coefficient')
#ax.set_title('')
ax.set_xticks(ind + width / 2)
sep = '~'
fontsize = 18;
ax.legend((rects1[0], rects2[0]), ('Repression', 'Activation'), frameon=False, fontsize=fontsize, loc ='upper left')
ax.set_xticklabels(('0h%s0h' % sep, '0h%s0.5h' % sep, '0.5h%s0.5h' % sep, '0h%s4h' % sep, '0.5h%s4h' % sep, '4h%s4h' % sep), rotation = 45)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(fontsize)

if(args.plot):
    plt.savefig(args.plot, format = os.path.basename(args.plot).split(".")[-1])
else:
    plt.show()




    
    
