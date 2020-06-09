#! /usr/bin/python
'''Draws replicated kinetics'''

import argparse
import sys
import os
from collections import defaultdict

#import pandas as pd;
import numpy as np;
from scipy.stats import sem
import matplotlib.pyplot as plt;

#from afbio.sequencetools import sliding_window
from afbio.numerictools import get_accumulated_derivatives,  smooth_with_averaging



parser = argparse.ArgumentParser(description='Draws replicated kinetics');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the kinetics table");
parser.add_argument('--concentration', nargs = '?', required=True, type = str, help = "Path to the concentration table");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory for the plots");
parser.add_argument('--format', nargs = '?', required=False, default = 'svg', type = str, help = "Format of the plots");
parser.add_argument('--ylimit', nargs = '?', const=True, default = False, type = bool, help = "If set, universal limit for Y-axis is used");
parser.add_argument('--derivative_window', nargs = '?', default=10, type = int, help = "Depth of the derivative");
parser.add_argument('--smooth_window', nargs = '?', default=10, type = int, help = "Flank length of the data smoothing");
parser.add_argument('--min_fraction', nargs = '?', default=0.05, type = float, help = "Min fraction of the max growth speed to set a threshold for start and end of an exponential phase");
parser.add_argument('--norm', nargs = '?', default=False, const=True, type = bool, help = "If set, data are normalized to the maximal intensity for each replicate individually");
args = parser.parse_args();

# set globals
COLORS = ['palevioletred', 'mediumvioletred', 'darksalmon']
FONTSIZE, LINEWIDTH = 28, 4

    

    
    

def merge_replicates(data):
    return data.mean(axis=0);


### PARSING FUNCTIONS ###

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




### GLOBAL PARAMETERS ###

def get_statistics(xvalues, yvalues, window, min_fraction):
    derivatives = get_accumulated_derivatives(xvalues, yvalues, window)
    growth_max_value = max(derivatives)
    limit = min_fraction*growth_max_value
    
    growth_max_index = np.argmax(derivatives);
    growth_start = [x[0] for x in enumerate(derivatives) if x[1]>=limit][0]
    growth_end = [x[0] for x in enumerate(derivatives[growth_max_index:]) if x[1]<limit][0] + growth_max_index
    
    data_start, data_end = growth_start+window, growth_end+window
    growth_total_speed = (yvalues[data_end]-yvalues[data_start])/(xvalues[data_end]-xvalues[data_start])
    plateau_value = np.mean(yvalues[data_end:data_end+window*5])
    acceleration_max = max(get_accumulated_derivatives(xvalues[window:-window], derivatives, window))
    
    
    return xvalues[growth_start+window], xvalues[growth_max_index + window], xvalues[growth_end + window], growth_max_value, growth_total_speed, plateau_value, acceleration_max


def process_sample_with_replicates(data, times, window, min_fraction):

    merged = merge_replicates(np.array(data))
    if( not any(merged)):
        return None
    
    local_list = [];
    for yvalues in data: 
        growth_start, growth_max_index, growth_end, growth_max_value, growth_total_speed, plateau_value, acceleration_max = get_statistics(times, yvalues, window, min_fraction)
        local_list.append(( growth_start, growth_max_index-growth_start, growth_end-growth_start, growth_total_speed, growth_max_value, plateau_value, (growth_max_index-growth_start)/(growth_end-growth_start), acceleration_max ))
    
    #print(len(local_list))
    local_list = np.array(local_list)
    return local_list.mean(axis=0), sem(local_list, axis=0), merged




### REPLICATES CONVERGENCY ###

def findstart(xvalues, yvalues, window, min_fraction):
    derivatives = get_accumulated_derivatives(xvalues, yvalues, window)
    limit = min_fraction*max(derivatives)
    return [x[0] for x in enumerate(derivatives) if x[1]>=limit][0]

#def smooth_with_averaging(vals, window):
    #derivatives = get_accumulated_derivatives(xvalues, yvalues, window)
    #limit = min_fraction*max(derivatives)
    #return [x[0] for x in enumerate(derivatives) if x[1]>=limit][0]

def converge_replicates(data, times, derivative_window, min_fraction, smooth_window):
    starts = [findstart(times, x, derivative_window, min_fraction) for x in data]
    aligned_data = [x[0][x[1]:] for x in zip(data, starts)]
    smoothed_data = [list(smooth_with_averaging(x, smooth_window)) for x in data]
    return smoothed_data, aligned_data
    


### GRAPHICS FUNCTIONS ###

def draw_technical(times, merged, sample, channel, format_, ylimit, data_limit=400):
    
    fig, ax1 = plt.subplots(figsize=(16,9))
    if(data_limit):
        times = times[:data_limit]
        merged = merged[:data_limit]
        
    ax1.plot(times, merged, color = COLORS[0], label = "signal") 
    ax1.spines['top'].set_visible(False)
    ax1.tick_params(axis='both', labelsize=FONTSIZE)
    plt.ylabel('ThT Fluorescence Intensity [a.u.]', fontsize = FONTSIZE)
    plt.xlabel('Time [h]', fontsize = FONTSIZE)    
    
    derivatives = get_accumulated_derivatives(times, merged, args.derivative_window)
    derivatives = [0]*args.derivative_window + derivatives + [0]*args.derivative_window
    ax2 = ax1.twinx()
    ax2.plot(times, derivatives, color = COLORS[1], label = "derivative")
    
    growth_start, growth_max_index, growth_end, growth_max_value, growth_total_speed, plateau_value, acceleration_max = get_statistics(times, merged, args.derivative_window, args.min_fraction)
    ax2.axvline(growth_start)
    ax2.axvline(growth_max_index)
    ax2.axvline(growth_end)
    
    if(ylimit):
        ax1.set_ylim(0, ylimit);
    plt.legend(frameon=False)
    plt.savefig(os.path.join(args.outdir, "technical_%s_channel%d.%s" % (sample, channel, format_)), format = format_)
    plt.close()
    plt.clf()
    
    
def draw_convergent(times, converged_data, sample, channel, format_, data_limit=400):
    fig, ax = plt.subplots(figsize=(16,9))    
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', labelsize=FONTSIZE)
    plt.ylabel('ThT Fluorescence Intensity [a.u.]', fontsize = FONTSIZE)
    plt.xlabel('Time [h]', fontsize = FONTSIZE)    
    
    for datum, color in zip(converged_data, COLORS):
        ax.plot(times[:data_limit], datum[:data_limit], color = color) 
    plt.savefig(os.path.join(args.outdir, "converged_%s_channel%d.%s" % (sample, channel, format_)), format = format_)
    plt.close()
    plt.clf()    


    
### PARSING ###

with open(args.path, 'r') as f:
    lines = f.readlines()
    wells = lines[6].strip().split(",")
    samples = lines[7].strip().split(",")
    position2sample = dict( [ x for x in enumerate(samples[2:]) ] )
    
    channel2lines = defaultdict(list);

    
    for line in lines[8:-1]:
        a = line.strip().split(",");
        channel2lines[a[0]].append(line);


sample2concentration = [];        
with open(args.concentration, 'r') as f:
    for l in f:
        a = l.strip().split("\t")
        sample2concentration.append(( "Sample X%s" % a[0], float(a[1]) ))
        
        
        
        
### EXECUTION ###

statistic_dict =  defaultdict(list)
error_dict = defaultdict(list)
names = []

DATA_LIMIT = 400
for channel, lines in enumerate(channel2lines.values(), start=1):
    times, sample2data = proces_lines(lines);
    if(args.ylimit):
        temp = [];
        for data in sample2data.values():
            temp.append(max([max(x) for x in data]));
        ylimit = max(temp)*1.01
        sys.stderr.write("Y-limit for the channel %d is set to %d\n" % (channel, ylimit))
    else:
        ylimit = False
        
    for sample, _ in sample2concentration:
        names.append(sample);
        data = sample2data[sample]
        if(args.norm):
            data = [np.array(x)-min(x) for x in data]
            data = [x/max(x) for x in data]
        stat, errors, merged = process_sample_with_replicates(data, times, args.derivative_window, args.min_fraction)
        statistic_dict[channel].append(stat)
        error_dict[channel].append(errors)
        
        converged_data, stub = converge_replicates(data, times, args.derivative_window, args.min_fraction, args.smooth_window)
        
        draw_convergent(times, converged_data, sample, channel, args.format, data_limit=400)
        draw_technical(times, merged, sample, channel, args.format, ylimit, data_limit=400)
        

    
        
    
    
### DRAW GLOBAL PARAMETERS ###

stypes = ['growth_start', 'growth_acceleration_duration', 'growth_length_duration', 'growth_average_speed', 'growth_max_speed', 'plateau', 'acceleration_fraction', 'acceleration_max']
ylabels = ['time [h]', 'time [h]', 'time [h]', 'growth speed [h^-1]', 'maximal growth speed [h^-1]', 'plateau intensity [a.u]', 'ratio', 'maximal acceleration [h^-2]']

 
#print(names)
for channel, statistic_list in statistic_dict.items():
    statistic_list = np.array(statistic_list).transpose();
    error_list = np.array(error_dict[channel]).transpose();
    for stype, ylabel, statistic, errors in zip(stypes, ylabels, statistic_list, error_list):
        #print(len(statistic))
        #print(statistic)
        fig, ax1 = plt.subplots(figsize=(16,9))
        fig.suptitle(stype, fontsize=FONTSIZE)
        plt.tight_layout(rect=[0.1, 0.15, 0.9, 0.9])
        x_range = range(len(statistic))

        
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.tick_params(axis='both', labelsize=FONTSIZE)
        ax1.set_ylabel(ylabel, fontsize = FONTSIZE)
        ax1.set_xlabel('TH relative molar concentration', fontsize = FONTSIZE)
        ax1.set_xticks(x_range, minor=False)
        ax1.set_xticklabels(["None"] + ["%1.3f" % x[1] for x in sample2concentration[1:]], fontsize=FONTSIZE/2, rotation=90)
    
    
        #print(errors)
        ax1.plot(x_range, statistic, color = COLORS[0], linewidth=LINEWIDTH)
        ax1.errorbar(x_range, statistic, yerr=errors, linestyle='', capsize=12, linewidth=LINEWIDTH/2, capthick=LINEWIDTH/2, color='black')

        
        plt.savefig(os.path.join(args.outdir, "%s_channel%d.%s" % (stype, channel, args.format)), format = args.format)
        plt.close()
        plt.clf()
    
    

    
    
        


### DEPRICATED ###

#def get_derivatives(xvalues, yvalues):
    #d_xvals = [(x[1] + x[0])/2 for x in zip(xvalues[1:], xvalues)]
    #d_yvals = [(x[0] - x[1])/(x[2]-x[3]) for x in zip(yvalues[1:], yvalues, xvalues[1:], xvalues)]
    #return np.array(d_xvals), np.array(d_yvals)

### OLD VERSION ###
#def get_accumulated_derivatives(xvalues, yvalues, window):
    
    #d_xvals = [(x[1] + x[0])/2 for x in zip(xvalues[window*2-1:], xvalues)]
    #d_yvals = np.zeros(len(d_xvals))
    #for w in range(window):
        #l_yvals = [(x[0] - x[1])/(x[2]-x[3]) for x in zip(yvalues[2*w+1:], yvalues, xvalues[2*w+1:], xvalues)]
        #flank = window-w-1;
        #if(flank):
            #l_yvals = np.array(l_yvals[flank:-flank]);
        #d_yvals+=l_yvals;
    #return d_xvals, d_yvals/window
    
#def get_ratio_derivative(xvalues, yvalues):
    #d_xvals = [(x[1] + x[0])/2 for x in zip(xvalues[1:], xvalues)]
    #d_yvals = [(x[0]/x[1])/(x[2]-x[3]) if x[1] else 0 for x in zip(yvalues[1:], yvalues, xvalues[1:], xvalues)]
    #return d_xvals, d_yvals
    
    
#def get_stable_derivatives(xvalues, yvalues):
    #d_xvals = xvalues[1:-1]
    #d_yvals = [(x[0] - x[1])/(x[2]-x[3]) for x in zip(yvalues[1:], yvalues, xvalues[1:], xvalues)]
    #d_yvals = [np.mean(x) for x in zip (d_yvals[1:], d_yvals)]
    #return np.array(d_xvals), np.array(d_yvals)
        
    
        
