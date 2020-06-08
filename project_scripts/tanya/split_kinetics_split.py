#! /usr/bin/python
'''Draws replicated kinetics'''

import argparse
import sys
import os
from collections import defaultdict

#import pandas as pd;
import numpy as np;
import matplotlib.pyplot as plt;



parser = argparse.ArgumentParser(description='Draws replicated kinetics');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the kinetics table");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory for the plots");
parser.add_argument('--format', nargs = '?', required=False, default = 'svg', type = str, help = "Format of the plots");
parser.add_argument('--ylimit', nargs = '?', const=True, default = False, type = bool, help = "If set, universal limit for Y-axis is used");
args = parser.parse_args();


def proces_lines(lines):
    position2values = defaultdict(list)
    times = [];  
    
    for line in lines:
        a = line.strip().split(",");
        times.append(int(a[1])/3600)
        for pos, value in enumerate(a[2:]):
            if(value):
                position2values[pos].append(float(value))
            else:
                position2values[pos].append(0);
    
    sample2data = defaultdict(list)
    for pos, values in position2values.items():
        sample2data[position2sample[pos]].append(values);

    return times, sample2data


with open(args.path, 'r') as f:
    lines = f.readlines()
    wells = lines[6].strip().split(",")
    samples = lines[7].strip().split(",")
    position2sample = dict( [ x for x in enumerate(samples[2:]) ] )
    
    channel2lines = defaultdict(list);

    
    for line in lines[8:-1]:
        a = line.strip().split(",");
        channel2lines[a[0]].append(line);
        
        
        
        
colors = ['palevioletred', 'mediumvioletred', 'darksalmon']    


for channel, lines in enumerate(channel2lines.values(), start=1):
    times, sample2data = proces_lines(lines);
    #print(times)
    if(args.ylimit):
        temp = [];
        for data in sample2data.values():
            temp.append(max([max(x) for x in data]));
        ylimit = max(temp)*1.01
        sys.stderr.write("Y-limit for the channel %d is set to %d\n" % (channel, ylimit))
        
        
    for sample, data in sample2data.items():
        plt.figure(1, figsize=(16,10), frameon=False)
        ax = plt.subplot(111)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', labelsize='x-large')
        plt.ylabel('ThT Fluorescence Intensity [a.u.]', fontsize = 'x-large')
        plt.xlabel('Time [h]', fontsize = 'x-large')
        plt.suptitle("channel_%d.%s" % (channel, sample), fontsize=50)

        for c, (values, color) in enumerate(zip(data, colors)):
            ax.plot(times, values, color = color, label = "replicate_%d" % (c+1))
        if(args.ylimit):
            ax.set_ylim(0, ylimit);
        ## manually choose the xlim ##
        #ax.set_xlim(0,12)
        plt.legend(frameon=False)
        plt.savefig(os.path.join(args.outdir, "channel_%d.%s.%s" % (channel, sample, args.format)), format = args.format)
        plt.clf()
    
    
        
    
    
        
