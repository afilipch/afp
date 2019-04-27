import argparse
import os
import sys
import numpy as np;
import matplotlib.pyplot as plt;
import os
from collections import defaultdict
from matplotlib.patches import Rectangle

from pybedtools import BedTool, Interval
import pandas as pd

parser = argparse.ArgumentParser(description='Represents peaks evolution in course of time-series experiments');
parser.add_argument('--regions', nargs = '?', required=True, type = str, help = "Path to the consensus regions");
parser.add_argument('--coverage', nargs = '+', required=True, type = str, help = "Path to the coverage files, the order of the given files must correspond to the time series");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the genomic annotation, plots will be annotated with the provided genomic features");
parser.add_argument('--custom', nargs = '?', default=False, const=True, type = bool, help = "If set the annotation genes are supposed to be already processed, if not they are supposed to be in NCBI gff3 format");

parser.add_argument('--flen', nargs = '?', default=40, type = int, help = "Length of the peak\'s flanks to be drawn");
parser.add_argument('--overlap', nargs = '?', default=0.5, type = float, help = "Minimal reciprocal overlap fraction required for the peaks to be considered as the same peak");
parser.add_argument('--normalize', nargs = '?', default=False, const=True, type = str, help = "If set, normalized coverage is output");
args = parser.parse_args();

def getcoverage(path, normalize):
    coverage = pd.read_csv(path, sep="\t" , names = ["chr", "postion", "coverage"]).coverage.values
    if(normalize):
        coverage = coverage/np.mean(coverage);
    return coverage


regions = BedTool(args.regions)
coveragesets = [getcoverage(x, args.normalize) for x in args.coverage]


def splitcoverage(region, coveragesets, flank):
    signal = [np.concatenate((np.zeros(flank), x[region.start:region.end], np.zeros(flank))) for x in coveragesets]
    control  = [np.concatenate((x[max(0, region.start-flank):region.start], np.zeros(region.end - region.start), x[region.end:region.end +flank])) for x in coveragesets]
    return list(zip(signal, control))

#print(len(regions))
signal_noise_sets = [splitcoverage(x, coveragesets, args.flen) for x in regions];


### Plot functions ###################################################################################################
def setticks(start, end):
    l = float(end - start);
    a = [5,10,20,50,100,200,500,1000, 2000]
    for el in a:
        if( 3 < l/el < 9):
            scale = el;
            break;
    else:
        return [];
    
    
    s = start//scale*scale
    if(s < start):
        s += scale;
        
    e = end//scale*scale
    if(e == end):
        e += scale;
    
        
    locs = [x for x in range(s, e+scale, scale)]
    labels = [str(x) for x in locs]
    return locs, labels


def draw_annotation(ax, locan, hlim):
    for c, an in enumerate(locan):
        rect = Rectangle( (an[0], 0.25+c), an[1] - an[0], 0.5, facecolor = 'green', edgecolor = 'green')
        ax.add_patch(rect)
        

ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])    
    
if(args.normalize):
    ylabel = 'normalized read coverage (num of means)'
else:
    ylabel = 'read coverage'


def getplot(signal_noise_local, size, start, end, path, fignum, locan, zscore):
    fig, axes = plt.subplots(nrows=size, ncols = 1, sharex=False, sharey=True, figsize = (16, 5*(size+1)), frameon=False)
    plt.suptitle("Top %s peak with z-score %s" % (ordinal(fignum), zscore), fontsize = 'xx-large')
    plt.tight_layout(rect=[0.1, 0.15, 0.95, 0.96], h_pad = 4)
    xticklocs, xticklabels = setticks(start, end)
    for ax, (signal, noise) in zip(axes, signal_noise_local):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()   
        ax.tick_params(axis='both', labelsize='x-large')
        ax.bar(range(start, end), height = noise, width=1, linewidth = 0, color='lightblue')
        ax.bar(range(start, end), height = signal, width=1, linewidth = 0, color=(242/256, 97/256, 68/256))
        ax.set_xlim(start, end)
        ax.set_xticks(xticklocs)
        ax.set_xticklabels(xticklabels)
    
    ansize = len(locan)
    box = axes[-1].get_position()
    x0 = box.x0
    x1 = box.x1-x0
    ylen = box.y1 - box.y0
    anax = fig.add_axes([x0, 0.15 - ylen/2 - ylen/4, x1, ylen/2])
    anax.spines['top'].set_visible(False)
    anax.spines['right'].set_visible(False)
    anax.spines['left'].set_visible(False) 
    anax.tick_params(axis='both', labelsize='x-large')
    anax.set_ylim(0, ansize)
    anax.set_xlim(start, end)
    anax.set_xticks(xticklocs)
    anax.set_xticklabels(xticklabels)
    anax.set_yticks([0.5 + x for x in range(ansize)])
    anax.set_yticklabels([ "%s(%s)" % (x[2], x[3]) for x in locan])
    draw_annotation(anax, locan, max([max(x[0]) for x in signal_noise_local]) )
    
    
    fig.text(0.5, 0.01, 'genomic position (nt)', ha='center', fontsize='xx-large')
    fig.text(0.005, 0.5, ylabel, va='center', rotation='vertical', fontsize='xx-large')
    plt.savefig(path, format='png')
    #plt.show()
    plt.close()
    
##########################################################################################################################


regions = BedTool([Interval(x.chrom, max(0, x.start - args.flen), x.end + args.flen, name=x.name, score=x.attrs['zscores'] , strand=x.strand) for x in regions])



annotation = BedTool(args.annotation)
if(not args.custom):
    annotation = [x for x in annotation if x in ['gene', 'pseudogene']]



rawannotated = regions.intersect(annotation, wo = True)
annotated = defaultdict(list)

def get_annotation(intersection, offset):
    return max(intersection.start, int(intersection[offset+3])), min(intersection.end, int(intersection[offset+4])), dict( [x.strip().split('=') for x in intersection[offset+8].split(";")])['Name'], intersection[offset+6]

for el in rawannotated:
    an = get_annotation(el, 6)
    annotated[el.name].append(an)

    


size = len(signal_noise_sets[0])
#regions = list(regions);
#regions.sort(key = lambda x: float(x.score), reverse = True)
for c, (region, snl) in enumerate(zip(regions, signal_noise_sets)):
    locan = annotated[region.name];
    locan.sort(key=lambda x: x[3]);
    #sys.stdout.write(str(region))
    os.path.join(args.outdir, "peak%d" % (c+1))
    getplot(snl, size, region.start, region.end, os.path.join(args.outdir, "peak%d.png" % (c+1)), c+1, locan, max([float(x) for x in region.score.split(",") if x != 'None']))
    if(not (c+1) % 10):
        sys.stderr.write("%d peaks are processed\n" % (c+1))
    #sys.exit()





            
