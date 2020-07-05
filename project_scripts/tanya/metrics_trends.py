#! /usr/bin/python
'''Draws replicated kinetics'''

import argparse
import sys
import os
from collections import defaultdict
from pathlib import Path

#import pandas as pd;
import numpy as np;
from scipy.stats import sem
import matplotlib.pyplot as plt;
#from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from itertools import  zip_longest
#from afbio.sequencetools import sliding_window
from afbio.numerictools import get_accumulated_derivatives,  smooth_with_averaging



parser = argparse.ArgumentParser(description='Draws replicated kinetics');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the kinetics table");
parser.add_argument('--concentration', nargs = '?', required=True, type = str, help = "Path to the concentration table");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory for the plots");
parser.add_argument('--format', nargs = '?', required=False, default = 'svg', type = str, help = "Format of the plots");
parser.add_argument('--ylimit', nargs = '?', const=True, default = False, type = bool, help = "If set, universal limit for Y-axis is used");
parser.add_argument('--derivative_window', nargs = '?', default=1, type = int, help = "Depth of the derivative");
parser.add_argument('--smooth_window', nargs = '?', default=1, type = int, help = "Flank length of the data smoothing");
parser.add_argument('--minfraction', nargs = '?', default=0.05, type = float, help = "Min fraction of the max growth speed to set a threshold for start and end of an exponential phase");
parser.add_argument('--outliers', nargs = '+', default = [], type = str, help = "Outlier tracks to be removed. Must be provided in format \'[channnel]![sample]![replicate num]\'")
parser.add_argument('--norm', nargs = '?', default=False, const=True, type = bool, help = "If set, data are normalized to the maximal intensity for each replicate individually");
args = parser.parse_args();

# set globals
COLORS = ['palevioletred', 'mediumvioletred', 'darksalmon']
RAINBOW_COLORS = ["darkmagenta", "midnightblue", "cornflowerblue", "mediumturquoise", "seagreen", "yellowgreen", "darkgoldenrod", "tomato", "saddlebrown", "darkorange", "gold", "khaki"]
FONTSIZE, LINEWIDTH = 28, 4

#print(args.outliers)
outliers = [tuple(x.split("!")) for x in args.outliers]

#Path(args.outdir).mkdir(parents=True, exist_ok=True)

dir_technical = os.path.join(args.outdir, 'technical'); 
Path(dir_technical).mkdir(parents=True, exist_ok=True)

dir_replicates_converged = os.path.join(args.outdir, 'replicates_converged'); 
Path(dir_replicates_converged).mkdir(parents=True, exist_ok=True)

dir_replicates_real = os.path.join(args.outdir, 'replicates_real'); 
Path(dir_replicates_real).mkdir(parents=True, exist_ok=True)

dir_converged = os.path.join(args.outdir, 'converged'); 
Path(dir_converged).mkdir(parents=True, exist_ok=True)

dir_metrics = os.path.join(args.outdir, 'metrics'); 
Path(dir_metrics).mkdir(parents=True, exist_ok=True)


    

    
    

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

def get_statistics(xvalues, yvalues, window, minfraction):
    derivatives = get_accumulated_derivatives(xvalues, yvalues, window)
    growth_max_value = max(derivatives)
    limit = minfraction*growth_max_value
    
    growth_max_index = np.argmax(derivatives);
    growth_start = [x[0] for x in enumerate(derivatives) if x[1]>=limit][0]
    growth_end = [x[0] for x in enumerate(derivatives[growth_max_index:]) if x[1]<limit][0] + growth_max_index
    
    data_start, data_end = growth_start+window, growth_end+window
    growth_total_speed = (yvalues[data_end]-yvalues[data_start])/(xvalues[data_end]-xvalues[data_start])
    plateau_value = np.mean(yvalues[data_end:data_end+window*5])
    acceleration_max = max(get_accumulated_derivatives(xvalues[window:-window], derivatives, window))
    
    
    return xvalues[growth_start+window], xvalues[growth_max_index + window], xvalues[growth_end + window], growth_max_value, growth_total_speed, plateau_value, acceleration_max


def process_sample_with_replicates(data, times, window, minfraction):

    merged = merge_replicates(np.array(data))
    if( not any(merged)):
        return None
    
    local_list = [];
    for yvalues in data: 
        growth_start, growth_max_index, growth_end, growth_max_value, growth_total_speed, plateau_value, acceleration_max = get_statistics(times, yvalues, window, minfraction)
        local_list.append(( growth_start, growth_max_index-growth_start, growth_end-growth_start, growth_total_speed, growth_max_value, plateau_value, (growth_max_index-growth_start)/(growth_end-growth_start), acceleration_max ))
    
    #print(len(local_list))
    local_list = np.array(local_list)
    return local_list.mean(axis=0), sem(local_list, axis=0), merged




### REPLICATES CONVERGENCY ###

def findstart(xvalues, yvalues, window, minfraction):
    derivatives = get_accumulated_derivatives(xvalues, yvalues, window)
    limit = minfraction*max(derivatives)
    return [x[0] for x in enumerate(derivatives) if x[1]>=limit][0]

#def smooth_with_averaging(vals, window):
    #derivatives = get_accumulated_derivatives(xvalues, yvalues, window)
    #limit = minfraction*max(derivatives)
    #return [x[0] for x in enumerate(derivatives) if x[1]>=limit][0]

def converge_replicates(data, times, derivative_window, minfraction, smooth_window):
    starts = [findstart(times, x, derivative_window, minfraction) for x in data]
    common_start = int(np.mean(starts));
    tails = [  x[0][ max(0, x[1] - common_start): x[1] ] for x in zip(data, starts) ]
    tails = [x[::-1] for x in tails]
    common_tail = [np.mean([x for x in y if x != None]) for y in zip_longest(*tails)]
    common_tail = common_tail[::-1]
    #for l in zip_longest(*tails):
        #print(np.mean([x for x in l if x != None]))
    
    aligned_data = [x[0][x[1]:] for x in zip(data, starts)]
    minlen = min([len(x) for x in aligned_data])
    aligned_data = [common_tail + list(x[:minlen]) for x in aligned_data]
    
    arr = np.array(aligned_data)
    #errors = sem(arr, axis = 0)
    #print(len(errors))
    #sys.exit()
    return aligned_data, np.mean(arr, axis=0), sem(arr, axis = 0)
    


### GRAPHICS FUNCTIONS ###

def draw_technical(times, merged, window, sample, channel, format_, ylimit, outdir, data_limit=400):
    
    fig, ax1 = plt.subplots(figsize=(16,9))
    derivatives = get_accumulated_derivatives(times, merged, window)
    derivatives = [0]*window + derivatives + [0]*window
    
    growth_start, growth_max_index, growth_end, growth_max_value, growth_total_speed, plateau_value, acceleration_max = get_statistics(times, merged, window, args.minfraction)
    if(data_limit):
        times = times[:data_limit]
        merged = merged[:data_limit]
        derivatives = derivatives[:data_limit]

    ax1.plot(times, merged, color = COLORS[0], label = "signal") 
    ax1.spines['top'].set_visible(False)
    ax1.tick_params(axis='both', labelsize=FONTSIZE)
    plt.ylabel('ThT Fluorescence Intensity [a.u.]', fontsize = FONTSIZE)
    plt.xlabel('Time [h]', fontsize = FONTSIZE)    
    
    ax2 = ax1.twinx()
    ax2.plot(times, derivatives, color = COLORS[1], label = "derivative")
    
    ax2.axvline(growth_start)
    ax2.axvline(growth_max_index)
    ax2.axvline(growth_end)
    
    if(ylimit):
        ax1.set_ylim(0, ylimit);
    plt.legend(frameon=False)
    plt.savefig(os.path.join(outdir, "technical_%s_channel%d.%s" % (sample, channel, format_)), format = format_)
    plt.close()
    plt.clf()
    
    
def draw_convergent_replicates(times, converged_data, converged_mean, converged_sem, sample, channel, format_, outdir, data_limit=400):
    fig, ax = plt.subplots(figsize=(16,9))    
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', labelsize=FONTSIZE)
    plt.ylabel('ThT Fluorescence Intensity [a.u.]', fontsize = FONTSIZE)
    plt.xlabel('Time [h]', fontsize = FONTSIZE)    
    
    for datum, color in zip(converged_data, COLORS):
        ax.plot(times[:data_limit], datum[:data_limit], color = color) 
    #draw averaged trend with error shadows
    ax.plot(times[:data_limit], converged_mean[:data_limit], 'k', color='#3F7F4C')
    ax.fill_between(times[:data_limit], converged_mean[:data_limit]-converged_sem[:data_limit], converged_mean[:data_limit]+converged_sem[:data_limit], alpha=0.5, edgecolor='#3F7F4C', facecolor='#7EFF99', linewidth=0)
    
    #print(outdir)
    plt.savefig(os.path.join(outdir, "converged_%s_channel%d.%s" % (sample, channel, format_)), format = format_)
    plt.close()
    plt.clf() 
    
    
def draw_real_replicates(times, data, sample, channel, format_, outdir, data_limit=400):
    fig, ax = plt.subplots(figsize=(16,9))    
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', labelsize=FONTSIZE)
    plt.ylabel('ThT Fluorescence Intensity [a.u.]', fontsize = FONTSIZE)
    plt.xlabel('Time [h]', fontsize = FONTSIZE)    
    
    for datum, color in zip(data, COLORS):
        ax.plot(times[:data_limit], datum[:data_limit], color = color) 
    
    #print(outdir)
    plt.savefig(os.path.join(outdir, "real_%s_channel%d.%s" % (sample, channel, format_)), format = format_)
    plt.close()
    plt.clf() 
    
    
    

def draw_all_convergent(times, converged_mean_list, converged_sem_list, labels, channel, format_, outdir, colors = RAINBOW_COLORS, data_limit=400):
    fig, ax = plt.subplots(figsize=(16,9))    
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', labelsize=FONTSIZE)
    plt.ylabel('ThT Fluorescence Intensity [a.u.]', fontsize = FONTSIZE)
    plt.xlabel('Time [h]', fontsize = FONTSIZE)    
    
    #dl = data_limit + max(start_list)
    x = times[:data_limit]
    
    handles = []
    for converged_mean, converged_sem, color in zip(converged_mean_list, converged_sem_list, colors):
        y = converged_mean[:data_limit]
        errors = converged_sem[:data_limit]       
        ax.plot(x, y, 'k', color=color, linewidth = LINEWIDTH/2)
        ax.fill_between(x, y-errors, y+errors, alpha=0.5, edgecolor=color, facecolor=color, linewidth=0)
        #handles.append(Circle( (0.5,0.5), radius=0.1, facecolor=color, edgecolor=color) )
        handles.append(Line2D([0], [0], marker='o', color="white", markerfacecolor=color, markersize=FONTSIZE*0.8) )
    #ax.axhline(0, color="white", linewidth = LINEWIDTH)
    plt.legend(handles=handles, labels=labels, loc='upper left', frameon = False, fontsize=FONTSIZE*0.8)
    #print(outdir)
    plt.savefig(os.path.join(outdir, "all_converged_channel%d.%s" % (channel, format_)), format = format_)
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
    
    converged_mean_list = [];
    converged_sem_list = [];
    labels = [];
    
    #start_list = [];
    
    if(args.ylimit):
        temp = [];
        for data in sample2data.values():
            temp.append(max([max(x) for x in data]));
        ylimit = max(temp)*1.01
        sys.stderr.write("Y-limit for the channel %d is set to %d\n" % (channel, ylimit))
    else:
        ylimit = False
        
    for sample, concentration in sample2concentration:
        names.append(sample);
        data = sample2data[sample]
        
        #print(channel, sample.split(" ")[1])
        if(outliers):
            #print(outliers)
            #print((channel, sample.split(" ")[1], str(1)))
            #print()
            data = [x[1] for x in enumerate(data) if (str(channel), sample.split(" ")[1], str(x[0]+1)) not in outliers]
        
        
        if(args.smooth_window > 1):
            #print(data[0][-10:])
            data = [list(smooth_with_averaging(x, args.smooth_window))[:-2*args.smooth_window] for x in data]
            #print(data[0][-10:])
        data_limit = min(int(len(data[0])*0.6), DATA_LIMIT)
        draw_real_replicates(times, data, sample, channel, args.format, dir_replicates_real, data_limit=data_limit)
        
        if(args.norm):
            data = [np.array(x)-min(x) for x in data]
            data = [x/max(x) for x in data]

        stat, errors, merged = process_sample_with_replicates(data, times, args.derivative_window, args.minfraction)
        statistic_dict[channel].append(stat)
        error_dict[channel].append(errors)
        
        converged_data, converged_mean, converged_sem = converge_replicates(data, times, args.derivative_window, args.minfraction, args.smooth_window)
        converged_mean_list.append(converged_mean);
        converged_sem_list.append(converged_sem)
        labels.append("c_%1.4f" % concentration)
        #start_list.append(start)
        
        data_limit = min(int(len(converged_data[0])*0.6), DATA_LIMIT)
        draw_convergent_replicates(times, converged_data, converged_mean, converged_sem, sample, channel, args.format, dir_replicates_converged, data_limit=data_limit)
        draw_technical(times, merged, args.derivative_window, sample, channel, args.format, ylimit, dir_technical, data_limit=data_limit)
        
    draw_all_convergent(times, converged_mean_list, converged_sem_list, labels, channel, args.format, dir_converged, data_limit=data_limit)
        

    
        
    
    
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

        
        plt.savefig(os.path.join(dir_metrics, "%s_channel%d.%s" % (stype, channel, args.format)), format = args.format)
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
        
    
        
