#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores high peak in prophage covereage'''

import argparse
import os
import sys
import math
from collections import defaultdict, Counter, namedtuple
from itertools import product, combinations
from glob import glob
import copy


from Bio import Entrez
import numpy as np;

from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;
from matplotlib.patches import Circle
import matplotlib.patches as mpatches



parser = argparse.ArgumentParser(description='Explores high peak in prophage covereage');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the prophage table, gff format");
#parser.add_argument('--bacteria', nargs = '?', required=True, type = str, help = "Path to the bacterial annotation");
parser.add_argument('--maxlength', nargs = '?', default=10e10, type = int, help = "Maximal length allowed for prophages");
parser.add_argument('--plot', nargs = '?', required=True, type = str, help = "Path to the folder for graphical output");
parser.add_argument('--format', nargs = '?', default='eps', choices = ['eps', 'svg', 'png', 'pdf'], type = str, help = "Format of the graphical output");
parser.add_argument('--btype', nargs = '?', default='Bacteria', type = str, help = "Type of bacteria to be analyzed");
parser.add_argument('--maxfraction', nargs = '?', default=0.05, type = float, help = "Maximal allowed prophage genome length (as a fraction of host genome length) to be reported as circle");
args = parser.parse_args();

#Prophage = namedtuple('Prophage', ['bacteria', 'start', 'end', 'num_of_genes', 'strand', 'integrase'])


### Reading the input

prophages = BedTool(args.path)

# filter length
prophages_filtered = [x for x in prophages if len(x) <= args.maxlength]
sys.stderr.write("%d prophages passed the length filtering (length <= %d) out of %d\n\n" % (len(prophages_filtered), args.maxlength, len(prophages)))
prophages = prophages_filtered;
prophages_filtered = [x for x in prophages if args.btype in x.attrs['host_family']]
sys.stderr.write("%d %s prophages are selected from total %d prophages\n\n" % (len(prophages_filtered), args.btype, len(prophages)))
prophages = prophages_filtered;




### Circular plot

scale = 1000;


#def prophage2interval:
    #return int(prophage.attrs['host_start']), int(prophage.attrs['host_end'])
    
intervals_intact = [];
intervals_disabled = [];
for prophage in prophages:
    #print(prophage.attrs)
    interval = int(prophage.attrs['host_start']), int(prophage.attrs['host_end'])
    if(prophage.attrs['integrase'] == 'True'):
        intervals_intact.append(interval);
    else:
        intervals_disabled.append(interval);

sys.stderr.write("%d number of intact prophages\n%d number of disabled prophages\n\n" % (len(intervals_intact), len(intervals_disabled)))            

size_threshold = scale*args.maxfraction    
intervals_intact = [ x for x in intervals_intact if x[1] - x[0] < size_threshold]
intervals_disabled = [ x for x in intervals_disabled if x[1] - x[0] < size_threshold]
 
### Plot with coverage  
 
def intervals2coverage(intervals, scale):
    coverage = np.zeros(scale+1);
    for s, e in intervals:
        coverage[s:e] += 1;
    return coverage;
        
coverage_intact = intervals2coverage(intervals_intact, scale);
coverage_disabled = intervals2coverage(intervals_disabled, scale);
#coverage_disabled[:200] = 0
#coverage_disabled[201:] = 0


#set up the plot
bottom = 10
max_height = 8/max(max(coverage_intact), max(coverage_disabled));
pos2grade = 2*np.pi/(scale+1);

fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(polar=True))
#plt.tight_layout(rect=[0.24, 0.24, 0.75, 0.75], h_pad = 4)
#ax = plt.subplot(111, polar=True)
ax.grid(False)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.spines['polar'].set_visible(False)

theta = [-x*pos2grade+np.pi/2 for x in range(scale+1)]
width = pos2grade

radii_intact = [x*max_height for x in coverage_intact]
bars = ax.bar(theta, radii_intact, width=width, bottom=bottom, color='lightblue', linewidth=0)
radii_disabled = [x*max_height for x in coverage_disabled]
bars = ax.bar(theta, radii_disabled, width=width, bottom=np.repeat(bottom, scale+1)-radii_disabled, color='coral', linewidth=0)


legend = plt.legend(['intact prophages (%d)' % len(intervals_intact), 'prophage remnants (%d)' % len(intervals_disabled)], loc=(-0.1, 0.75), frameon=False, fontsize = 'xx-large')
ax.annotate("origin of replication", xy=(np.pi/2, bottom), xytext=(np.pi/2, bottom*1.5), arrowprops=dict(arrowstyle="->"), ha='center', fontsize = 'xx-large')
ax.plot(np.linspace(0, 2*np.pi, 100), [bottom]*100, color='r', linewidth=1)
plt.savefig(os.path.join(args.plot, '%s_%1.2f_prophage_coverage.%s' % (args.btype, args.maxfraction, args.format)), format=args.format)



### Plot with coverage  with circles


def collapse_intervals(intervals, maxshift):
    temp = copy.copy(intervals);
    collapsed = [];
    used = set();
    for i1 in intervals:
        #print()
        if(i1 not in used):
            length = i1[1] - i1[0];
            for i2 in temp:
                if( abs(i1[0] - i2[0])/length < maxshift and abs(i1[1] - i2[1])/length < maxshift ):
                    #print(i1,i2)
                    used.add(i2);
                    collapsed.append(i1);
    return Counter(collapsed);
        
    

def polar2cartesian(theta, rad):
    y = math.sin(theta)*rad
    x = math.cos(theta)*rad
    return x, y



def draw_circle(interval, scale, ax, occupied, inner, color, alpha):
    #print(interval)
    pos = (interval[1] + interval[0])/2
    size = ((interval[1] - interval[0])/scale)*np.pi*bottom
    theta = np.pi/2-pos2grade*pos
    if(inner):
        rad = bottom-size-occupied[int(pos)];
        occupied[interval[0]: interval[1]] += size
    else:
        rad = bottom+size+occupied[int(pos)];
        occupied[interval[0]: interval[1]] += size
    patch = plt.Circle(polar2cartesian(theta, rad), radius=size, color = color, transform=ax.transData._b, alpha=min(0.2+alpha/4, 1))
    ax.add_line(patch);
    #ax.autoscale()
    



#fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(polar=True))
#plt.tight_layout(rect=[0.24, 0.24, 0.75, 0.75], h_pad = 4)
ax = plt.subplot(111, polar=True)
ax.grid(False)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.spines['polar'].set_visible(False)


occupied = np.zeros(scale+1);
counter_intact = collapse_intervals(intervals_intact[:], 0.05)
print(len(counter_intact));
for interval, alpha in counter_intact.items():
    if(interval[1] - interval[0] < size_threshold):
        draw_circle(interval, scale, ax, occupied, False, 'lightblue', alpha)
        
        
occupied = np.zeros(scale+1);
counter_disabled = collapse_intervals(intervals_disabled[:], 0.05)
print(len(counter_disabled));
for interval, alpha in counter_disabled.items():
    draw_circle(interval, scale, ax, occupied, True, 'coral', alpha)


legend_patches = mpatches.Patch(color='lightblue', label='intact prophages (%d)' % len(intervals_intact)), mpatches.Patch(color='coral', label='prophage remnants (%d)' % len(intervals_disabled))
legend = plt.legend(handles=legend_patches, loc=(-0.1, 0.75), frameon=False, fontsize = 'xx-large')
ax.annotate("origin of replication", xy=(np.pi/2, bottom), xytext=(np.pi/2, bottom*1.5), arrowprops=dict(arrowstyle="->"), ha='center', fontsize = 'xx-large')
ax.annotate("terminus", xy=(3*np.pi/2, bottom), xytext=(3*np.pi/2, bottom*1.5), arrowprops=dict(arrowstyle="->"), ha='center', fontsize = 'xx-large')

ax.plot(np.linspace(0, 2*np.pi, scale+1), [bottom]*(scale+1), color='red', linewidth=1)
ax.plot(np.linspace(0, 2*np.pi, 100), [bottom*2]*100, color='coral', linewidth=0)
#draw_circle(intervals_intact[10], ax)
#ax.autoscale()
plt.savefig(os.path.join(args.plot, '%s_%1.2f_prophage_circles.%s' % (args.btype, args.maxfraction, args.format)), format=args.format)



    
#sys.exit()
### Length plot    
lengths_intact = [len(x) for x in prophages if x.attrs['integrase'] == 'True']
lengths_disabled = [len(x) for x in prophages if x.attrs['integrase'] == 'False']
threshold = 200000
brange = np.linspace(0,threshold, 41)
bars_intact = [len(list(filter(lambda y: x[0]<y<=x[1], lengths_intact))) for x in zip(brange, brange[1:])] + [len([x for x in lengths_intact if x>threshold])]
bars_disabled = [len(list(filter(lambda y: x[0]<y<=x[1], lengths_disabled))) for x in zip(brange, brange[1:])] + [len([x for x in lengths_disabled if x>threshold])]

#print(bars)
#print(brange)
fig, ax = plt.subplots(figsize=(16, 9))
plt.tight_layout(rect=[0.04, 0.1, 0.95, 0.96], h_pad = 4)
ax.bar(np.arange(len(brange)), bars_intact, 0.25, color='lightblue')
ax.bar(np.arange(len(brange))+0.25, bars_disabled, 0.25, color='coral')

brange = np.linspace(0,threshold, 11)

plt.xticks(np.linspace(0,40, 11), [int(x) for x in brange[:-1]] + [">%d" % threshold], fontsize = 'xx-large', rotation=30)

plt.ylabel("number of prophages", fontsize = 'xx-large')
plt.xlabel("size [bp]", fontsize = 'xx-large')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize('xx-large')
legend = plt.legend(['intact prophages (%d)' % len(lengths_intact), 'prophage remnants (%d)' % len(lengths_disabled)], loc=(0.6, 0.8), frameon=False, fontsize = 'xx-large')
plt.savefig(os.path.join(args.plot, '%s_prophage_length.%s' % (args.btype, args.format)), format=args.format)





### Num of phages per bacteria plot    
brange = np.arange(1,12)
prophages_per_bacteria_intact = Counter([x.chrom for x in prophages if x.attrs['integrase'] == 'True']).values()
prophages_per_bacteria_disabled = Counter([x.chrom for x in prophages if x.attrs['integrase'] == 'False']).values()

def get_bars(prophages_per_bacteria):
    bars = np.zeros(11)
    for v in prophages_per_bacteria:
        if(v>10):
            bars[10] += 1;
        else:
            bars[v-1] += 1;
    return bars;

bars_intact = get_bars(prophages_per_bacteria_intact)
bars_disabled = get_bars(prophages_per_bacteria_disabled)
            
    
fig, ax = plt.subplots(figsize=(16, 9))
plt.tight_layout(rect=[0.04, 0.1, 0.95, 0.96], h_pad = 4)
ax.bar(brange, bars_intact, 0.25, color='lightblue')
ax.bar(brange+0.25, bars_disabled, 0.25, color='coral')
#plt.title(title)
plt.xticks(brange, [str(x) for x in brange[:-1]] + [">10"], fontsize = 'xx-large')
plt.ylabel("number of bacterial species", fontsize = 'xx-large')
plt.xlabel("number of prophages per single specie", fontsize = 'xx-large')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize('xx-large')
#plt.show()
legend = plt.legend(['intact prophages (%d)' % sum(prophages_per_bacteria_intact), 'prophage remnants (%d)' % sum(prophages_per_bacteria_disabled)], loc=(0.6, 0.8), frameon=False, fontsize = 'xx-large')
plt.savefig(os.path.join(args.plot, '%s_prophage_number_per_bacteria.%s' % (args.btype, args.format)), format=args.format)



















### Check for integrases
#viables = [];
#for c, fpath in enumerate(glob(args.path)):
    #viable = False;
    #with open(fpath) as f:
        #for l in f:
            #a =  l.strip().split();
            #fann = a[4];
            #if('integrase' in fann):
                #viable = True;
                #break;
    #viables.append(viable)
    #if(c and c % 1000 == 0):
        #print("%d files processed" % c);
#print(Counter(viables))