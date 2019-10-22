#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Annotates the discovered peaks'''

import argparse
import sys
import os
from collections import defaultdict
from bisect import bisect_right, bisect_left


import pandas as pd;
import numpy as np;
import matplotlib.pyplot as plt;
from pybedtools import BedTool

from afbio.pybedtools_af import construct_gff_interval
from afbio.sequencetools import coverage2dict

parser = argparse.ArgumentParser(description='Annotates the discovered peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the detected peaks");
parser.add_argument('--coverage', nargs = '?', required=True, type = str, help = "Path to the coverage track, bed format");
parser.add_argument('--convolution', nargs = '?', required=True, type = str, help = "adfsa");
args = parser.parse_args();




if(os.stat(args.path).st_size != 0):
    coverage_dict = coverage2dict(args.coverage);
    convolution = list(coverage2dict(args.convolution).values())[0]
    for peak in BedTool(args.path):
        top = int(peak.name)
        topcoverage = coverage_dict[peak.chrom][top]
        total_coverage = coverage_dict[peak.chrom][peak.start:peak.stop]
        newint = construct_gff_interval(peak.chrom, peak.start, peak.stop, 'binding_peak', score=peak.score, strand=peak.strand, source='un', frame='.', attrs=[("Name", peak.name), ("topcoverage", topcoverage), ("total_coverage", total_coverage)] )
        if(max(total_coverage) >= 1.1*topcoverage):
            fontsize = 24
            fig, ax1 = plt.subplots(figsize=(16,9))
            plt.tight_layout(rect=[0.1, 0.1, 0.95, 0.95])
            

            ax1.plot(total_coverage, 'b-')
            ax1.plot(top - peak.start, topcoverage, 'r*', linewidth=5)
            ax1.set_xlabel("position (nt)", fontsize=fontsize)
            ax1.set_ylabel('coverage', color='b', fontsize=fontsize)
            ax1.tick_params('y', colors='b')
            
            ax2 = ax1.twinx()
            ax2.plot(convolution[peak.start:peak.stop], 'g-')
            ax2.set_ylabel("convolution", color='g', fontsize=fontsize)
            ax2.tick_params('y', colors='g')
            
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False) 
            ax1.tick_params(axis='both', which='major', labelsize=fontsize)
            ax1.tick_params(axis='both', which='minor', labelsize=fontsize)
            plt.show()

        
        
    
    


    







