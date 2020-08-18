#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Draws a heatmap for glxr expression'''

import argparse
import os
import sys
from scipy.stats import pearsonr
from itertools import combinations
import numpy as np;
import pandas as pd;
import matplotlib.pyplot as plt;
from math import log
from collections import defaultdict


parser = argparse.ArgumentParser(description='Draws a heatmap for glxr expression');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the expression table");
parser.add_argument('--plot', nargs = '?', default='', type = str, help = "Path for the output plot");
parser.add_argument('--top', nargs = '?', default=50, type = int, help = "Number of top genes to take");
args = parser.parse_args();

ORDER = ['glu_wt', 'glu_ko_cyab', 'ace_glu_wt', 'ace_glu_ko_cyab']

def merge_with_replicates(exp_labels, mylist):
    res = defaultdict(list)
    for label, val in zip(exp_labels, mylist):
        key = "_".join(label.split('_')[:-1])
        res[key].append(float(val));
    res = [log(np.mean(res[x])+1, 2) for x in ORDER]
    #res = [(x[0], log(x[1]+1, 2)) for x in res.items()]
    #res.sort(key=lambda x: x[0])
    return res
    


### Read the input
expr_start = 4;
gene_pos = 2;
gene_expression = []
#gene_labels = []
with open(args.path) as f:
    exp_labels = [x.split(" ")[1] for x in next(f).strip().split("\t")[expr_start:]]
    #replicates_dict = dict([ (x[0], "_".join(x[1].split('_')[:-1])) for x in enumerate(exp_labels) ])
    #print(replicates_dict)
    for l in f:
        a = l.strip().split("\t")
        #gene_labels.append(a[gene_pos]);
        #expression.append([ log(float(x)+1, 2) for x in a[expr_start:]])
        gene_expression.append( (a[gene_pos], merge_with_replicates(exp_labels, a[expr_start:])) )

gene_expression.sort(key = lambda x: sum(x[1]), reverse = True)
gene_expression = gene_expression[:args.top]

expression = [x[1] for x in gene_expression]
gene_labels = [x[0] for x in gene_expression]
cmatrix = np.array(expression)



######################################################################################
### Draw a heatmap

#cmap="GnBu"
cmap="inferno_r"
fontsize = 22
fig, ax = plt.subplots(figsize=(8, 30))
plt.tight_layout(rect=[0.05, 0, 0.9, 0.9])
im = ax.imshow(cmatrix, cmap=cmap)
cbar = ax.figure.colorbar(im, ax=ax, cmap=cmap)
cbar.ax.set_ylabel('Log2(Peak Intensity + 1)', rotation=-90, va="bottom", fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

ax.set_xticks(np.arange(len(ORDER)))
ax.set_yticks(np.arange(len(gene_labels)))
ax.set_xticklabels(ORDER, rotation=90, fontsize=fontsize)
ax.set_yticklabels(gene_labels, fontsize=fontsize)
ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

for edge, spine in ax.spines.items():
    spine.set_visible(False)

ax.set_xticks(np.arange(cmatrix.shape[1]+1)-.5, minor=True)
ax.set_yticks(np.arange(cmatrix.shape[0]+1)-.5, minor=True)
#ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
ax.tick_params(which="minor", bottom=False, left=False)

if(args.plot):
    _format = args.plot.split(".")[-1]
    plt.savefig(args.plot, format = _format)
else:
    plt.show()
    
    
#
#I developed and coded a vastly novel algorithm to detect physiologically relevant molecular interactions. The traces of such interactions were hidden in publicly available datasets for years. Me and my colleague were the first to discover them. These traces were rare and weak, which set high demands on the algorithm I had to develop for their detection. I managed to balance on a sharp edge between sensitivity and specificity, designing novel approaches for not yet encountered problems.  

#For the last two years I have been working on 4-6 independent projects at the same time. I would not claim any of these projects to be extremely challenging. But manage 5 projects at the same time is a huge challenge itself. I advanced a lot in an art of setting right priorities and time management.  

#Tech skills: I deeply studied math in my college. Then switched to physics in the University. My PhD I got in the field of molecular biology and bioinformatics. I have more than 10 years experience in programming and more than 5 in data science. Therefore I believe that I can handle the problems at the crossroads of multiple fields. I believe that there should be a great idea which requires expertise in math, biology, physics and data science, and I will be a perfect fit to translate it into reality. 

#Business skills: My farther is running his own small company and I was raised in an atmosphere of entrepreneurship. I was expected to step for family business and got some practical skills. However I selected a path of science.

#It was at my PhD interview. I was applying to one of the most prominent system biology labs in Europe (and competitive to the best institutions in US). There were 150 other applicants to this single position. I decided to make something unusual and invested many efforts into visual and entertaining aspects of my presentation. I do not know, whether this strategy helped me, but I was selected.

#If you want to make something well you should learn from the best ones. I want to make a startup and I hope EF will give an opportunity to learn from the most successful startup entrepreneurs.

 #I am also would be happy to meet talented people with a various background in the cohort. I think such an exchange will help to develop new ideas, which I would never develop on my own.
