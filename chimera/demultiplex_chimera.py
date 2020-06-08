#! /usr/lib/python
'''Assignes to each read, if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits'''
import argparse
import os;
import sys;

import pysam;


from afbio.chimera import arwlist2chimeras
from afbio.generators import generator_segments




parser = argparse.ArgumentParser(description='Assignes to each read if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
#parser.add_argument('-o', '--output', nargs = '?', default = "sam", type = str, help = "path to the output folder");
#parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for output files, should reflect nature of mapping reference");

parser.add_argument('-mg', '--maxgap', nargs = '?', default = 4, type = int, help = "Max gap allowed between parts of chimera. It is also used to calculate chimera score. The more the maxgap, the less important is the gap.");
parser.add_argument('-mo', '--maxoverlap', nargs = '?', default = 6, type = int, help = "Max overlap allowed between parts of chimera. It is also used to calculate chimera score. The more the maxoverlap, the less important is the overlap.");
parser.add_argument('-sd', '--splice_distance', nargs = '?', default = 10000, type = int, help = "If the distance between chimeric parts on a reference is less than a splice_distance, then this chimera will be preferentially selected as potential splice junction");

parser.add_argument('--s_distance', nargs = '?', default = 12, type = float, help = "minimal distance allowed between the best and the second best hit. If the actual distance is less, than hit will be assigned as nonunique");
parser.add_argument('--ch_distance', nargs = '?', default = 12, type = float, help = "minimal distance allowed between the best and the second best chimera. If the actual distance is less, than chimera will be assigned as nonunique");

parser.add_argument('--reassign', nargs = '?', default = False, const=True, type = bool, help = "If set, coordinates will be reassigned");
args = parser.parse_args();


def select_splicing(chimeras, splice_distance):
    return [x for x in chimeras if (x.coordinates[0][0] == x.coordinates[1][0] and x.coordinates[0][1] == x.coordinates[1][1] and abs(x.coordinates[0][2] - x.coordinates[1][2])<=splice_distance)]


def chimera2score(chimera, gapbase, overlapbase):
    '''Score function for Chimera based on alignment score and gap'''
    if(chimera.gap>0):
        gap = chimera.gap
        overlap = 0
    else:
        gap = 0;
        overlap = -chimera.gap

    return chimera.AS*(1 - gap/gapbase)*(1-overlap/overlapbase)


def get_chimeras(arwlist, maxgap, maxoverlap, chimera_distance, splice_distance, reassign):
    chimeras = arwlist2chimeras(arwlist, maxgap, maxoverlap)
    if(not chimeras):
        return None

    for chimera in chimeras:
        chimera.get_coordinates(reassign)
        chimera.set_score(chimera2score(chimera, maxgap, maxoverlap))
        
    bestscore = max([x.score for x in chimeras])
    chimeras = [x for x in chimeras if x.score+chimera_distance>bestscore]

    spliced = select_splicing(chimeras, splice_distance)
    if(spliced):
        chimeras = spliced;


    if(len(chimeras) == 1):
        return chimeras[0];
    else:
        return None;


def get_single(arws, chimera, s_distance):
    if(chimera):
        chscore = chimera.score
    else:
        chscore = 0;
        
    bestscore = max([x.score for x in arws])

    if(chscore > bestscore+s_distance):
        return chimera, None
    else:
        narws = [x for x in arws if x.score+s_distance>bestscore]
        if(len(narws)==1):
            return None, narws[0];
        else:
            return None, None;

#mappings derived from the same read are pulled together. Collapsed into one (or more, in case of non-unique mappings) wrapping read object. 
#then they are written to a new destination, according to their source: real, or control


#with open(os.path.join(args.output, "%s.chimera.bed" % args.name), "w") as f:
for c, arwlist in enumerate(generator_segments(args.path)):
    chimera = get_chimeras(arwlist, args.maxgap, args.maxoverlap, args.ch_distance, args.splice_distance, args.reassign)
    chimera, single = get_single(arwlist, chimera, args.s_distance)
    if(chimera):
        sys.stdout.write('%s\n' % chimera.doublebed());
            
        #if(c and c % 100000 == 0):
            #sys.stderr.write('Reads processed:\t%d\n' % c);


#sys.stderr.write('\nChimeras detected:\t%d\nChimeras control detected:\t%d\nSingles detected:\t%d\nSingles control detected:\t%d\n' % tuple(counts));


