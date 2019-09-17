'''here we discriminate binding peaks from the other genome'''
import argparse
import os
import sys
import copy
import numpy as np;
from collections import defaultdict, Counter
from pybedtools import BedTool, Interval
from Bio import SeqIO
import matplotlib.pyplot as plt;

from afbio.sequencetools import get_at_content, sliding_window, coverage2dict
from afbio.pybedtools_af import construct_gff_interval

#from afbio.sequencetools import get_at_content, sliding_window, array2fixed_length, coverage2dict




parser = argparse.ArgumentParser(description='here we discriminate binding peaks from the other genome');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks, gff format");
parser.add_argument('--tracks', nargs = '+', required=True, type = str, help = "Path to AT content tracks");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts coordinates, gff file");
parser.add_argument('--pflank', nargs = '?', default=60, type = int, help = "Length of the flank around a peak center");
parser.add_argument('--thresholds', nargs = '+', required=True, type = float, help = "AT content thresholds, the order must correspond to the --tracks");

parser.add_argument('--clear', nargs = '?', default=100, type = int, help = "Additional flank length around the peak (after adding pflank) to remove from the no-peak regions");
parser.add_argument('--long_flank', nargs = '?', default=1000, type = int, help = "Flank length to calculate AT around the peaks");

parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");

parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
args = parser.parse_args();

def num_intersected_peaks(center_dict, areas, distance):
    count = 0
    local_dict = defaultdict(list)
    for area in areas:
        local_dict[area.chrom].append((area.start+area.stop)//2)
    for chrom, cl in center_dict.items():
        ll = local_dict[chrom]
        for c in cl:
            d = min([abs(c-x) for x in ll]) 
            if(d<= distance):
                count+=1
    return count;


def get_at_flanks(interval, genome, flank, margin):
    f1 = get_at_content(genome[interval.chrom][interval.start-flank:interval.start-margin].seq)
    f2 = get_at_content(genome[interval.chrom][interval.stop+margin:interval.stop+flank].seq)
    return (f1+f2)/2


def get_upstreams(transcripts, ulen, dlen):
    upstreams = []
    
    for transcript in transcripts:
        if(transcript.strand == '-'):
            start = transcript.stop - dlen
            stop = transcript.stop + ulen;
        else:
            start = transcript.start - ulen
            stop = transcript.start + dlen;
        if(start > 0):
            upstreams.append(Interval(transcript.chrom, start, stop, transcript.name, transcript.score, transcript.strand))
    return BedTool(upstreams)


def mask_track(track_dict, length):
    for arr in track_dict.values():
        arr[:length] = 0
        arr[-length:] = 0;


def get_closest(pos, positions):
    a = [abs(pos-x) for x in positions];
    i = np.argmin(a)
    return positions[i], pos - positions[i]

def area2interval(area, pflank, clear):
    chrom = area[0][4]
    start = area[0][0]
    stop = area[-1][0]+1
    score = str(max([x[3] for x in area]))
    distance = int(abs((np.mean([x[2] for x in area]))))
    peak = area[0][1]
    if(distance<pflank):
        _type = 'in'
    elif(distance<clear):
        _type = 'unk'
    else:
        _type = 'out'
    
    return construct_gff_interval(chrom, start, stop, 'long_at', score=score, strand='.', source='un', frame='.', attrs=[("peak", peak), ("distance", distance), ("type", _type), ("ID", str((stop+start)//2))])
    


def find_closest_peak(track_dict, center_dict, threshold):
    areas = [];
    for chrom, track in track_dict.items():
        area = [];
        centers = center_dict[chrom]
        curpos = -1;
        for pos, at in enumerate(track):
            if(at > threshold):
                center, distance = get_closest(pos, centers)
                if(pos - curpos > 40):
                    if(area):
                        areas.append(area);
                    area = [(pos, center, distance, at, chrom),]
                else:
                    area.append((pos, center, distance, at, chrom))
                curpos = pos;
        else:
            area.append((pos, center, distance, at, chrom))
                #print(pos, center, distance, at)
    return areas;



def find_closest_area(area, areas):
    temp = [x for x in areas if x.chrom==area.chrom and x.start!=area.start]
    distances = [abs(area.start-x.start) for x in temp];
    i = np.argmin(distances)
    return temp[i], distances[i]
    
    




genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
transcripts = BedTool(args.transcripts);
upstreams = get_upstreams(transcripts, 140, 20)
regions = BedTool(args.path)
center_dict = defaultdict(list);
for region in regions: 
    center_dict[region.chrom].append( (region.end+region.start)//2 )
#regions = [convert_region(x, args.pflank) for x in regions]
#regions = BedTool([x for x in regions if float(x.score)>args.zscore])
#regions_dict = dict([ (x.name, x) for x in regions ])


long_track_dict = coverage2dict(args.tracks[0])
mask_track(long_track_dict, 200)

areas = find_closest_peak(long_track_dict, center_dict, args.thresholds[0])
areas = BedTool([area2interval(x, args.pflank, args.clear) for x in areas])
temp_areas = []
for area in areas.intersect(b=upstreams, f=0.5, u=True):
    area.attrs['upstream'] = 'True'
    temp_areas.append(area);
for area in areas.intersect(b=upstreams, f=0.5, v=True):
    area.attrs['upstream'] = 'False'
    temp_areas.append(area)

areas = [];
for area in temp_areas:
    c_area, distance = find_closest_area(area, temp_areas)
    area.attrs['at_flank'] = str(get_at_flanks(area, genome, args.long_flank, 200))
    area.attrs['c_area'] = c_area.name
    area.attrs['c_area_distance'] = str(distance);
    areas.append(area)

areas.sort(key=lambda x: float(x.score));
areas = BedTool(areas)

for at_flank_thresh in np.linspace(0.42, 0.50, 9):
    areas = BedTool([x for x in areas if float(x.attrs['at_flank'])>=at_flank_thresh])
    
    print(at_flank_thresh)
    print()
    areas_in = BedTool([x for x in areas if x.attrs['type']=='in'])
    areas_out = BedTool([x for x in areas if x.attrs['type']=='out'])
    areas_unk = BedTool([x for x in areas if x.attrs['type']=='unk'])
    print("\nTotal peaks: %d\nPeaks overlapping long AT areas: %d\n" % (len(regions), num_intersected_peaks(center_dict, areas, args.pflank) ) )
    print("Total AT long areas: %d\nAT areas inside peaks: %d\nAT areas away from peaks: %d\nAT areas close to peaks: %d\n" % tuple([len(x) for x in (areas, areas_in, areas_out, areas_unk)]))
    print("#"*140)







sys.exit()
#############################################################################################################
#############################################################################################################





### UPSTREAM REGIONS PART
def get_upstream_fraction(areas):
    n = len([x for x in areas if x.attrs['upstream']=='True'])
    l = len(areas)
    if(l):
        return n/l*100
    else:
        return 0


print("Inside peaks AT areas upstream fraction: %1.1f%%\nAway from peak AT areas upstream fraction: %1.1f%%\nClose to peaks AT areas upstream fraction: %1.1f%%\n" %  tuple([get_upstream_fraction(x) for x in (areas_in, areas_out, areas_unk)]))


distance_thresholds = [100, 200, 500, 1000, 2000, 5000, 10000]
('Distance\tInside peaks\tClose to peaks\tAway from peaks')
for dt in distance_thresholds:
    a = [dt];
    for t_areas in (areas_in, areas_out, areas_unk):
        a.append(len([x for x in t_areas if int(x.attrs['c_area_distance'])<dt])/len(t_areas)*100)
    print('<%d\t%1.1f%%\t%1.1f%%\t%1.1f%%' % tuple(a))
print()    
    
    
distance_thresholds = [0, 100, 200, 500, 1000, 2000, 5000, 10000, 1000000]
('Distance Range\tFraction of Areas\tUpstream of them')
for dt1, dt2 in zip(distance_thresholds[:-1], distance_thresholds[1:]):
    passed_areas = [x for x in areas_out if int(x.attrs['c_area_distance'])>=dt1 and  int(x.attrs['c_area_distance'])<dt2 ]
    a = dt1, dt2, len(passed_areas), get_upstream_fraction(passed_areas)

    print('%d~%d\t%d\t%1.1f%%' % a)
    
    
#percentiles = [10, 20, 30, 40, 50, 60, 70, 80, 90]
#for t_areas in (areas_in, areas_out, areas_unk):
    #distances = [int(x.attrs['c_area_distance']) for x in t_areas]
    #for p in percentiles:
        #d_threshold = np.percentile(distances, p) - 0.0001;
        #print(p, d_threshold, len([x for x in distances if x>d_threshold]))
    #print()
    
    
    
    
#for area in areas_out:
    #if(int(area.attrs['c_area_distance']) > 2000):
        #print("%1.2f\t%s" % (float(area.score), area.attrs['upstream']) )
        
    

#for area in areas_unk:
    #c_area, distance = find_closest_area(area, areas)
    







    
    









    
