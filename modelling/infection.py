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
parser.add_argument('--moi', nargs = '?', required = True, type = int, help = "Multiplicity of infection");
parser.add_argument('--length_cycle', nargs = '?', default = 1, type = int, help = "Time of the cycle for the modelling in minutes, default: 1");
parser.add_argument('--max_cycles', nargs = '?', default = 200, type = int, help = "Number of cycles to model, default: 200");
#parser.add_argument('--cds', nargs = '?', type = str, help = "Path to the cds regions, gff format");
args = parser.parse_args();

def infection_basic(cell_number, cell_resistant, phage_number, cell_dead, infection_efficiency, parameters):
    cell_free = cell_number - cell_resistant
    cell_infected = int(cell_free*phage_number*infection_efficiency)
    cell_lysogenic = int(cell_infected*parameters['lysogenic'])
    cell_lysed = cell_infected - cell_lysogenic
    phage_produced = int(cell_lysed*parameters['phage_production'])

    return cell_number - cell_lysed, cell_resistant + cell_lysogenic, cell_dead + cell_lysed, phage_number + phage_produced, cell_lysogenic
    
    
def doubling_basic(cell_number, cell_resistant, food_limit, doubling_factor):
    return int(cell_number*doubling_factor), int(cell_resistant*doubling_factor), int(food_limit - cell_number*(doubling_factor-1))
    
 
#list_cell_number = [];
#list_cell_resistant = [];
#list_phage_number = [];
#list_food_limit = [];
def cycle_basic(cell_number, cell_resistant, phage_number, food_limit, cell_dead, doubling_factor, infection_efficiency, max_cycles, parameters):
    data = [(cell_number, cell_resistant, phage_number, cell_dead, food_limit),]
    for time in range(max_cycles):
        if(cell_number>cell_resistant):
            cell_free, cell_resistant, cell_dead, phage_number, cell_lysogenic = infection_basic(cell_number, cell_resistant, phage_number, cell_dead, infection_efficiency, parameters)
        if(food_limit>cell_number*(doubling_factor-1)):
            cell_number, cell_resistant, food_limit = doubling_basic(cell_number, cell_resistant, food_limit, doubling_factor)
            #print(food_limit)
        #if(time>320 and time < 360):
            #print(cell_resistant, cell_lysogenic)
        data.append((cell_number, cell_resistant, phage_number, cell_dead, food_limit))
    return np.array(data)
        
    


#Pre-set parameters
parameters = {"cell_number" : 10**6,
"doubling_time" : 40,
"infection_efficiency": 10**(-9),
"lysogenic" : 0.0001,
"defective_lysis" : 0.2,
"cell_resistant": 0.0001,
"lag_bio" : 20,
"phage_production": 5,
"food_limit" : 10**9}



#Set free parameters
infection_efficiency = parameters["infection_efficiency"]*args.length_cycle
doubling_factor = np.e**(np.log(2)/(parameters['doubling_time']/args.length_cycle))
print(doubling_factor)
                           


#Setting initials
food_limit = parameters['food_limit']
cell_number = parameters['cell_number']
cell_resistant = cell_number*parameters['cell_resistant']
phage_number = args.moi*cell_number
cell_dead = 0;

 
data = cycle_basic(cell_number, cell_resistant, phage_number, food_limit, cell_dead, doubling_factor, infection_efficiency, args.max_cycles, parameters)

 
 
################################ Drawing Section ################################ 

fontsize = 18;
labels = ["num of phages", "num of dead cells", "food availible"]
fig, axes = plt.subplots(nrows=2, ncols = 2, figsize = (16, 9), frameon=False)
(ax1, ax2), (ax3, ax4) = axes

ax1.set_xlabel('Time [minutes]', fontsize=fontsize)
ax1.set_ylabel("num of total cells", fontsize=fontsize)    
ax1.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
ax1.spines['top'].set_visible(False)
ax1.plot(data[:,0])
ax1.plot(data[:,1], color = 'orange')

#ax2 = ax1.twinx()
#ax2.set_ylabel("num of resistant cells", fontsize=fontsize)    
#ax2.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
#ax2.spines['top'].set_visible(False)
#ax2.plot(data[:,1])



axes = (ax2, ax3, ax4)
for count, (ax, label) in enumerate(zip(axes, labels), start=2):
    ax.set_xlabel('Time [minutes]', fontsize=fontsize)
    ax.set_ylabel(label, fontsize=fontsize)    
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.plot(data[:,count])
    
plt.show()

    
    
    
    
    
    






