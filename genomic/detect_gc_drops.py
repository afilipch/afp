#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Detects  GC drops along the provided GC content track'''

import argparse
import os
import sys
#import scipy
import numpy as np;
import pandas as pd;

from afbio.filters import dsk, usk, trk;
from afbio.peaks import convolute


parser = argparse.ArgumentParser(description='Detects  GC drops along the provided GC content track');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the GC content file, bed3 format");
#parser.add_argument('--kernel', nargs = '?', choices = ['ndg', 'square'], default='ndg', type = str, help = "Kernel to use for the fitering");
#parser.add_argument('--peakwidth', nargs = '?', default=0, type = int, help = "Expected width of the peaks, if not set peakwidth is estimated automatically");
#parser.add_argument('--widthfactor', nargs = '?', default=1, type = float, help = "Multiplier for the eastimated peakwidth, so peakwidth=estimated_peakwidth*widthfactor");
parser.add_argument('--lookup', nargs = '?', default=10, type = int, help = "Size of an area for finding local extrema in convolution track");
parser.add_argument('--bandwidth', nargs = '?', default=60, type = int, help = "Assumed mininmal length of a dip");
parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");
#parser.add_argument('--convolution', nargs = '?', default='', type = str, help = "If set, the convolution track is written to the provided file");
args = parser.parse_args();

###Read Coverage
gc_content = pd.read_csv(args.path, sep="\t" , names = ["chr", "position", "gc"]).gc.values


dkernel = dsk(up=1, down=-1)
ukernel = usk(up=1, down=-1)
#kernel = trk(steps=1001);
coverage = gc_content[:]#[1]*100 + [1/x for x in range(1, 101)] + [1]*100
#print(ukernel)

def find_local_extrema(convolution, lookup):
    extrema = [];
    for pos, c in enumerate(convolution[lookup: -lookup], start = lookup):
        if(c >= max(convolution[pos-lookup:pos]) and c > max(convolution[pos+1:pos + lookup+1]) ):
            extrema.append((pos, c, 'max'));
        elif(c <= min(convolution[pos-lookup:pos]) and c < min(convolution[pos+1:pos + lookup+1]) ):
            extrema.append((pos, c, 'min'));
    return extrema;


def find_local_minima(convolution, lookup):
    minima = [];
    for pos, c in enumerate(convolution[lookup: -lookup], start = lookup):
        if(c <= min(convolution[pos-lookup:pos]) and c < min(convolution[pos+1:pos + lookup+1]) ):
            minima.append((pos, c));
    return minima;



#sys.stderr.write("bb\n")
convolution_down = np.array(convolute(coverage, dkernel, args.bandwidth, threads=8, takepeak=False))
extrema = find_local_extrema(convolution_down, args.lookup);

for e in extrema:
    print("chr\t%d\t%1.4f\t%s" % e);


#print(maxima)

#print(coverage)
#convolution_up = np.array(convolute(coverage, ukernel, bandwidth, threads=8, takepeak=False))



###Plot coverage vs convolution
if(args.plot):
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(coverage, 'b-')
    ax1.set_xlabel("position (nt)")
    ax1.set_ylabel('coverage', color='b')
    ax1.tick_params('y', colors='b')
    
    #convolution_down = [max(0,x) for x in convolution_down]
    minima = [x[:2] for x in extrema if x[2] == 'min']
    ax2 = ax1.twinx()
    ax2.plot(convolution_down, 'r-')
    ax2.plot(*zip(*minima), "g*")
    ax2.set_ylabel("convolution", color='r')
    ax2.tick_params('y', colors='r')
    ax2.axhline(color='g')
    
    #convolution_up = [max(0,x) for x in convolution_up]
    #ax3 = fig.add_subplot(313)
    #ax3.plot(convolution_up, 'g-')
    #ax3.set_ylabel("convolution", color='g')
    #ax3.tick_params('y', colors='g')
    
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False) 

    fig.tight_layout()
    
    _format = args.plot.split(".")[-1]
    plt.show();
    #plt.savefig(args.plot, format = _format)




#limit = gc_content.mean()*0.9
#sys.stderr.write("limit is set to %1.5f\n\n" % limit);

#drops = [];

#indrop = False
#for pos, gcc in enumerate(gc_content):
    #if(gcc < limit and not indrop):
        #indrop = True;
        #start = pos;
    #elif(gcc >= limit and indrop):
        #indrop = False
        #end = pos;
        #drops.append((start, end));
#else:
    #if(indrop):
        #drops.append((start, len(gc_content)));
        
#for s, e in drops[:20]:
    #print("%d\t%d\t%s" % (s, e, gc_content[s-1: e+1]));
    
    
#print(np.percentile([x[1]-x[0] for x in drops], range(101)))
#print(max([x[1]-x[0] for x in drops]))
