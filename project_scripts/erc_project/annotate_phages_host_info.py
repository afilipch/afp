#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Adds additional host annotation to phages. Converts phages table from tsv to gff format'''

import argparse
import os
import sys
from glob import glob
from collections import defaultdict, namedtuple

from afbio.pybedtools_af import construct_gff_interval;





parser = argparse.ArgumentParser(description='Adds additional host annotation to phages. Converts phages table from tsv to gff format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the prophage table, tsv format");
parser.add_argument('--bacteria', nargs = '?', required=True, type = str, help = "Path to the bacterial annotation");
args = parser.parse_args();

Prophage = namedtuple('Prophage', ['bacteria', 'start', 'end', 'num_of_genes', 'strand', 'integrase'])


### Reading the input

prophages = [];
with open(args.path) as f:
    for l in f:
        a = l.strip().split("\t")
        a[1:4] = [int(x) for x in a[1:4]]
        if(a[5] == 'True'):
            a[5] = True;
        else:
            a[5] = False
        prophages.append(Prophage(*a));


bacterial_annotation = [];
with open(args.bacteria) as f:
    for l in f:
        a = l.strip().split("\t")
        a[1] = "".join(a[1].split());
        a[1] = ",".join(a[1].split(';'))
        a[2] = int(a[2])
        bacterial_annotation.append(a);


### Circular plot

bacteria2prophages = defaultdict(list);

for prophage in prophages:
    bacteria2prophages[prophage.bacteria].append(prophage)




for bname, bannotation, blength in bacterial_annotation:
    for prophage in sorted(bacteria2prophages[bname], key = lambda x: x.start):
        scale = 1000/blength
        host_start = round(prophage.start*scale)
        host_end = round(prophage.end*scale)+1;
        interval = construct_gff_interval(bname, prophage.start, prophage.end, 'prophage', score='0', strand=prophage.strand, source='ovidiu', frame='.', attrs=[('host_length', str(blength)),  
        ('host_family', bannotation), ('Name', "_".join([str(x) for x in (bname, prophage.start, prophage.end)])), ('rlength', "%1.4f" % ((prophage.end-prophage.start)/blength) ), 
        ('integrase', str(prophage.integrase)), ('host_start', "%d" % host_start ),  ('host_end', "%d" % host_end)])
        if(host_start>0 and host_end<1000):
            sys.stdout.write(str(interval));
        else:
            sys.stderr.write("\nProphage region is outside bacterial borders\n%s" % interval);
            
            
        
        