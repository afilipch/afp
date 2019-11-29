#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores how the region of interest changes the chap coverage in varieing conditions'''

import argparse
import os
import sys
import numpy as np;
import pandas as pd;
from collections import defaultdict, Counter
from pybedtools import BedTool
from afbio.sequencetools import coverage2dict

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


parser = argparse.ArgumentParser(description='Explores how the region of interest changes the chap coverage in varieing conditions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the directory with annotated binding peaks");
parser.add_argument('--interest', nargs = '?', required=True, type = str, help = "Path to the bed file with regions of interest, bed file");
parser.add_argument('--coverage', nargs = '?', required=True, type = str, help = "Path to the folder with coverage track, bed format")
#parser.add_argument('--smooth', nargs = '?', default=0, type = int, help = "Sliding window half-length used to smooth the control genomic coverage, default: 0");
#parser.add_argument('--plot', nargs = '?', type = str, help = "Path to the plot");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output plot directory")
args = parser.parse_args();

############################# PEAKS SECTION #############################

badship = 'glu_his_rep1';
goodship = 'glu_his_rep2';

def fix_data(data):
    for p, (name, value) in enumerate(data):
        if(name == badship):
            b = p;
        if(name == goodship):
            g = p;
    data[b] = data[g] 

peaks_list = [(x.split(".")[0],  BedTool(os.path.join(args.path, x)) ) for x in os.listdir(args.path) if 'annotated' in x];
interest = BedTool(args.interest);

multidata = [];
for pos in range(4):
    data = []
    for name, peaks in peaks_list:
        region = max(peaks.intersect(b=interest, u=True, f=0.35, F=0.35), key = lambda x: float(x.score))
        data.append(( name, float(region.attrs['other_coverage'].split(",")[pos]) ));
    data.sort(key= lambda x: x[0])
    fix_data(data)
    data = [  ( "_".join(x[0][0].split("_")[:2]), x[0][1], x[1][1]) for x in zip(data[::2], data[1::2])]
    data.sort(key= lambda x: x[0][::-1])


    multidata.append(data)
    
    
    
############################# COVERAGE SECTION #############################   
boxes = [(307, 580, 'NCgl1683'), (696, 710, 'GntR BS'), (766, 795, 'NCgl1682')]



flank = 400;
interval = interest[0];
x_range = np.arange(interval.start - flank, interval.stop + flank)
#cov_list = [x for x in os.listdir(args.path) if 'normalized' in x];
data = []
t = 0;
for x in os.listdir(args.coverage):
    if('normalized' in x):
        #t+=1
        name = tuple(x.split(".")[0].split("_"))
        cov_dict = coverage2dict(os.path.join(args.coverage, x), cpos = 4)
        coverage = cov_dict[interval.chrom][interval.start - flank : interval.stop + flank]
        #coverage = np.arange(10)*(t+1)
        data.append((name, coverage))
        
data.sort(key=lambda x: "".join(x[0]))
        
gntr_glu_data = [x for x in data if x[0][0] == 'glu' and x[0][1] == 'his']
gntr_glu_data = gntr_glu_data[1:]
gntr_fru_data = [x for x in data if x[0][0] == 'fru' and x[0][1] == 'his']
cgps_glu_data = [x for x in data if x[0][0] == 'glu' and x[0][1] == 'strep']
cgps_fru_data = [x for x in data if x[0][0] == 'fru' and x[0][1] == 'strep']

        
        
        
#cov_dict = coverage2dict(args.coverage, cpos = 4)
#interval = interest[0];
#coverage = cov_dict[interval.chrom][interval.start - flank : interval.stop + flank]

    
    
############################# DRAWING SECTION #############################
def draw(data, ctype, width = 0.6, fontsize = 24, linewidth = 4):
    bars = [(x[1]+x[2])/2 for x in data]
    errors = np.array([(x[1] - min(x[0][1:]), max(x[0][1:]) - x[1]) for x in zip(data, bars)]).transpose()


    fig, ax = plt.subplots(figsize = (16, 9))
    plt.tight_layout(rect=[0.08, 0.08, 1, 1])
    x = np.arange(len(data))
    ax.bar(x, bars, width, label='bound genes', color = 'lightblue')
    plt.errorbar(x, bars, yerr=errors, linestyle='', capsize=12, linewidth=linewidth/2, capthick=linewidth/2, color='black')
    #ax.bar(x+width, [x[2] for x in data], width, label='bound genes', color = 'lightblue')

    # add some text for labels, title and axes ticks
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('Name of experiment', fontsize=fontsize)
    ax.set_ylabel('Peak Intensity [averaged coverage]', fontsize=fontsize)
    ax.set_xticks(x)
    ax.set_xticklabels([x[0] for x in data], rotation = 0)
    ax.tick_params(axis='both', labelsize=fontsize-4, top=False, right=False)

    plt.savefig(os.path.join(args.outdir, "interesting_intensity_%s.%s") % (ctype, args.format) , format = args.format)
    
    
def draw_coverage(protein, glu, fru, x_range, boxes, fontsize = 24, linewidth = 4):
    fig, ax = plt.subplots(figsize = (16, 9))
    plt.tight_layout(rect=[0.08, 0.08, 1, 0.95])
    plt.title(protein, fontsize=fontsize+4)
    
    colors = 'lightblue', 'coral'
    markers = '-', '--'
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('position', fontsize=fontsize)
    ax.set_ylabel('Normalized coverage', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    
    ylimit = 0
    for dataset, color in zip([glu, fru], colors):
        for data, marker in zip(dataset, markers):
            ax.plot(x_range, data[1], color=color, linewidth=linewidth, linestyle=marker, label="%s_%s" % (data[0][0], data[0][2]))
            ylimit = max(ylimit, max(data[1]))
            
    for box in boxes:
        ax.add_patch(Rectangle( (box[0], 0), box[1]-box[0] , ylimit/18, color = 'lightgray'))
        if(box[2] == "GntR BS"):
            ax.text(box[0]-70 , 0 , box[2], color = 'black', fontsize=fontsize-4)
        else:
            ax.text(box[0] , 0 , box[2], color = 'black', fontsize=fontsize-4) 
            
    fig.legend(frameon=False, loc='upper right', fontsize=fontsize)
    plt.savefig(os.path.join(args.outdir, "interesting_coverage_%s.%s") % (protein, args.format) , format = args.format)
    
    

        
        
############################# EXECUTING SECTION #############################
for ctype, data in zip(['all', 'signal', 'noise', 'raw'], multidata):
    draw(data, ctype)
    
draw_coverage('gntr', gntr_glu_data, gntr_fru_data, x_range, boxes, fontsize = 24, linewidth = 4)
draw_coverage('cgps', cgps_glu_data, cgps_fru_data, x_range, boxes, fontsize = 24, linewidth = 4)
    
    

    


