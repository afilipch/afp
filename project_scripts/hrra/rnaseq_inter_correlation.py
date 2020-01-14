#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates table of the correlation between rna-seq replicates for each time point'''

import argparse
import os
import sys
import numpy as np;
from scipy.stats import pearsonr




parser = argparse.ArgumentParser(description='Draws a scatter plot peak TSS distances versus corresponding genes expression change');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the rnaseq differential tables");
#parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Output destination for the plot");
args = parser.parse_args();

print("Time-point\tWT\tKO")
for path in args.path:
    with open(path) as f:
        next(f);
        wt_list = []
        ko_list = []
        for l in f:
            a = l.strip().split("\t")
            wt = [float(x) for x in a[1].split(";")] 
            ko = [float(x) for x in a[2].split(";")]
            if(sum(wt)>5):
                wt_list.append(wt)
            if(sum(ko)>5):
                ko_list.append(ko)
            
        ko_list = np.array(ko_list);
        wt_list = np.array(wt_list);
        ko_corr = pearsonr(ko_list[:,1], ko_list[:,0])[0]
        wt_corr = pearsonr(wt_list[:,1], wt_list[:,0])[0]
        
        name = os.path.basename(path).split(".")[0]
        print("%s\t%1.4f\t%1.4f" % (name, wt_corr, ko_corr))
        
        
#print(wt_list[:,1])

