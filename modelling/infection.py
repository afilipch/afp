#! /usr/bin/python
'''Explores relationships between infection and growth'''

import argparse
import sys
import os
from collections import defaultdict

import numpy as np;
import matplotlib.pyplot as plt;



parser = argparse.ArgumentParser(description='Explores relationships between infection and growth');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Growth data from biolector");
parser.add_argument('--moi', nargs = '+', required = True, type = float, help = "Multiplicity of infection");
parser.add_argument('--length_cycle', nargs = '?', default = 1, type = int, help = "Time of the cycle for the modelling in minutes, default: 1");
parser.add_argument('--max_cycles', nargs = '?', default = 200, type = int, help = "Number of cycles to model, default: 200");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
#parser.add_argument('--cds', nargs = '?', type = str, help = "Path to the cds regions, gff format");
args = parser.parse_args();

def infection_basic(cell_number, cell_resistant, phage_number, cell_dead, food_limit, infection_efficiency, parameters):
    cell_free = cell_number - cell_resistant
    cell_infected = min(cell_free, int(cell_free*phage_number*infection_efficiency))
    cell_lysogenic = int(cell_infected*parameters['lysogenic'])
    cell_lysed = cell_infected - cell_lysogenic
    phage_produced = int(cell_lysed*parameters['phage_production'])

    return cell_number - cell_lysed, cell_resistant + cell_lysogenic, cell_dead + cell_lysed, phage_number + phage_produced, cell_lysogenic, food_limit - int(cell_lysed*parameters['lysis_food'])
    
    
def doubling_basic(cell_number, cell_resistant, food_limit, doubling_factor):
    return int(cell_number*doubling_factor), int(cell_resistant*doubling_factor), int(food_limit - cell_number*(doubling_factor-1))
    
 
#list_cell_number = [];
#list_cell_resistant = [];
#list_phage_number = [];
#list_food_limit = [];
def cycle_basic(cell_number, cell_resistant, phage_number, food_limit, cell_dead, doubling_factor, infection_efficiency, max_cycles, parameters):
    data = [(cell_number, cell_resistant, phage_number, cell_dead, food_limit),]
    for time in range(max_cycles):
        if(food_limit>cell_number*(doubling_factor-1)):
            if(cell_number>cell_resistant):
                cell_free, cell_resistant, cell_dead, phage_number, cell_lysogenic, food_limit = infection_basic(cell_number, cell_resistant, phage_number, cell_dead, food_limit, infection_efficiency, parameters)
            cell_number, cell_resistant, food_limit = doubling_basic(cell_number, cell_resistant, food_limit, doubling_factor)
            #print(food_limit)
        #if(time>320 and time < 360):
            #print(cell_resistant, cell_lysogenic)
        data.append((cell_number, cell_resistant, phage_number, cell_dead, food_limit))
    return np.array(data)
        
    


#Pre-set parameters
parameters = {"cell_number" : 10**6,
"doubling_time" : 40,
"infection_efficiency": 1*(10**(-10)),
"lysogenic" : 0.001,
"defective_lysis" : 0.2,
"cell_resistant": 0.001,
"lag_bio" : 20,
"phage_production": 5,
"lysis_food": 1,
"food_limit" : 0.5*(10**10)}



#Set free parameters
doubling_factor = np.e**(np.log(2)/(parameters['doubling_time']/args.length_cycle))
#print(doubling_factor)
                           




 
 
################################ Drawing Section ################################ 

def draw_info(data, moi, outdir, fontsize=18):
    labels = ["num of phages", "num of dead cells", "food availible"]
    fig, axes = plt.subplots(nrows=2, ncols = 2, figsize = (16, 9), frameon=False)
    (ax1, ax2), (ax3, ax4) = axes

    ax1.set_xlabel('Time [minutes]', fontsize=fontsize)
    ax1.set_ylabel("num of total cells", fontsize=fontsize)    
    ax1.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax1.spines['top'].set_visible(False)
    ax1.plot(data[:,0])
    ax1.plot(data[:,1], color = 'orange')



    axes = (ax2, ax3, ax4)
    for count, (ax, label) in enumerate(zip(axes, labels), start=2):
        ax.set_xlabel('Time [minutes]', fontsize=fontsize)
        ax.set_ylabel(label, fontsize=fontsize)    
        ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.plot(data[:,count])
        
    plt.savefig(os.path.join(outdir, "info.moi_%1.1f_.%s"  %  (moi, args.format) ) , format = args.format)
    plt.clf()
    
    
def draw_moi(total_list, outdir, ie_readable, fontsize=18):
    colors = plt.cm.RdPu(np.linspace(0.25,1,len(total_list)))
    #print(colors)
    fig, ax = plt.subplots(figsize=(16,9))

    ax.set_xlabel('Time [minutes]', fontsize=fontsize)
    ax.set_ylabel("num of total cells", fontsize=fontsize)    
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for total, moi, color in zip(total_list, args.moi, colors):
        ax.plot(total, color = color, label= "MOI=%1.1f" % moi)
    fig.legend(loc=(0.15, 0.75), frameon=False, fontsize=fontsize)   
        
    plt.savefig(os.path.join(outdir, "varied_moi.parameter_%s.%s"  %  (ie_readable, args.format) ) , format = args.format)
    plt.clf()
    plt.close()
    
def draw_parameter_distance(plateau_distances, ie_range, outdir, xlabel, fontsize=18):
    #colors = plt.cm.RdPu(np.linspace(0.25,1,len(plateau_distances)))
    fig, ax = plt.subplots(figsize=(16,9))

    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel("Discrimination between MOI", fontsize=fontsize)    
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    x_range = [np.log2(x) for x in ie_range]
    ax.plot(x_range, plateau_distances, color = 'darkblue')
    ax.plot(x_range, plateau_distances, 'r*')
    #fig.legend(loc=(0.15, 0.75), frameon=False, fontsize=fontsize)   
        
    plt.savefig(os.path.join(outdir, "parameter_discrimination.%s"  %  args.format ) , format = args.format)
    plt.clf()
    plt.close()
    
    
    
    
################################ Main Section ################################
def run(parameters, outdir, ie_readable='constant', if_info = True):
    total_list = [];
    for moi in args.moi:
        infection_efficiency = parameters["infection_efficiency"]*args.length_cycle
        food_limit = parameters['food_limit']
        cell_number = parameters['cell_number']
        cell_resistant = cell_number*parameters['cell_resistant']
        cell_dead = 0;
        phage_number = moi*cell_number
        #print(phage_number)
        data = cycle_basic(cell_number, cell_resistant, phage_number, food_limit, cell_dead, doubling_factor, infection_efficiency, args.max_cycles, parameters)
        total_list.append(data[:,0])
        if(if_info):
            draw_info(data, moi, outdir, fontsize=18)
    
    draw_moi(total_list, outdir, ie_readable, fontsize=18)
    return total_list;


run(parameters, args.outdir);


################################ Vary infection infection efficiency ################################
def get_plateau_diff(total_list):
    return [max(x[0]) - max(x[1]) for x in zip(total_list, total_list[1:])]

plateau_distances = []
ie_start = 10**(-11)
ie_range = np.array([2**x for x in range(11)])
#ie_range = np.array([1, 2, 5, 10, 100, 1000, 10000]);
for ie_cur in ie_range:
    ie_local = ie_cur*ie_start
    parameters['infection_efficiency'] = ie_local;
    total_list = run(parameters, os.path.join(args.outdir, 'efficiency'), ie_readable=ie_cur,  if_info=False);
    plateau_distances.append(min([abs(x) for x in get_plateau_diff(total_list)]));
    
draw_parameter_distance(plateau_distances, ie_range, os.path.join(args.outdir, 'efficiency'), 'Log2(Infection efficiency)', fontsize=18)
parameters['infection_efficiency'] = ie_start*2;

################################ Vary food availibility ################################
plateau_distances = []
food_start = 10**7
food_range = np.array([2**x for x in range(21)])
for f_cur in food_range:
    f_local = f_cur*food_start
    parameters['food_limit'] = f_local
    total_list = run(parameters, os.path.join(args.outdir, 'food'), ie_readable=f_cur,  if_info=False);
    plateau_distances.append(min([abs(x) for x in get_plateau_diff(total_list)]));
    
draw_parameter_distance(plateau_distances, food_range, os.path.join(args.outdir, 'food'), 'Log2(Food Limit)', fontsize=18)
    





#print( "\t".join(["TIME/MOI"] + [str(x) for x in args.moi]) )
#curtime = 0;
#for cycle in zip(*total_list):
    #print(cycle)
    

    
    
    
    
    
    






