#! /usr/bin/python
'''Collapses double bed/gff file of chimeric reads into interactions. That is, script merges chimeras with simultaneously intersecting regions'''

import sys
import argparse;
from collections import defaultdict

from pybedtools import BedTool

from afbio.pybedtools_af import construct_gff_interval

parser = argparse.ArgumentParser(description='converts bed-like file of chimeric reads into interactions. That is merging chimeras with intersecting regions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to chimeras, double bed/gff format. NOTE: Intervals are assumed to be sorted with sort.py, which is crucial");
parser.add_argument('-d', '--distance', nargs = '?', default = -1, type = int, help = "Minimum overlap(negative number)/maximum distance(positive number) of any pair of intervals in nucleotides to merge them");
parser.add_argument('-n', '--name', nargs = '?', default = 'i', type = str, help = "Base name for interactions")
parser.add_argument('-od', '--dictionary', nargs = '?', required = True, type = str, help = "Path to the output \"interaction id to read id\" file")
args = parser.parse_args();


def merge(intervals, name):
    chrom = intervals[0].chrom
    start = min([int(x.start) for x in intervals])
    stop = max([int(x.stop) for x in intervals])
    score = str(max([float(x.score) for x in intervals]))
    strand = intervals[0].strand

    gap = str(min([int(x.attrs['gap']) for x in intervals], key = abs))
    unique_reads = sum([int(x.attrs.get('n_uniq', 1)) for x in intervals])
    chscore = str(max([x.attrs['chscore'] for x in intervals]))

    return construct_gff_interval(chrom, start, stop, 'interaction', score=score, strand=strand, source='chiflex', frame='.', attrs=[("ID", name), ('gap', gap), ('chscore', chscore), ('n_uniq', unique_reads)])




def generate_overlaping_intervals(intervals, overlap):
    '''From a given iterable of intervals yields lists of intervals intersecting each other with each pairwise overlap more then a given cutoff. NOTE: The function is applicable only for sorted(with strandness) iterbales of intervals'''	
    first = intervals[0];
    merged = [first]
    start, stop = first.start, first.stop
    ref = (first.chrom, first.strand)

    for interval in intervals[1:]:
        cref = (interval.chrom, interval.strand)
        if(ref == cref):
            if(min(interval.stop, stop) - max(interval.start, start) >= overlap):
                merged.append(interval);
                stop = max(interval.stop, stop)
            else:
                yield merged;
                start, stop = interval.start, interval.stop;
                merged = [interval];
        else:
            yield merged
            ref = cref;
            start, stop = interval.start, interval.stop
            merged = [interval];
    yield merged


bedfile = BedTool(args.path)

if(len(bedfile) == 0):
    sys.stderr.write("input file is empty\n");
    sys.exit();


cid2iid = defaultdict(list);
for iid, merged_intervals in enumerate(generate_overlaping_intervals(bedfile, -1*args.distance)):
    for mi in merged_intervals:
        cid2iid[mi.name.split("|")[0]].append(iid);

            
iid2intervals = defaultdict(list);
for interval in bedfile:
    iid = tuple(sorted(cid2iid[interval.name.split("|")[0]]));
    iid2intervals[iid].append(interval);


with open(args.dictionary, 'w') as odf:
    for c, intervals in enumerate(iid2intervals.values()):
        
        name = "%s_%d" % (args.name, c)
        name_left = "%s|0" % name
        name_right = "%s|1" % name
        
        intervals_left = [x for x in intervals if x.name.split("|")[-1] == '0']
        sys.stdout.write(str(merge(intervals_left, name_left)));
        
        intervals_right = [x for x in intervals if x.name.split("|")[-1] == '1']
        sys.stdout.write(str(merge(intervals_right, name_right)));
        
        for li in intervals_left:		
            odf.write("%s\t%s\n" % (li.name.split("|")[0], name))
        
        
        
