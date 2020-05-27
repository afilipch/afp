#! /usr/lib/python
'''produce chimeras from merged and already filtered sam file'''
import argparse;
import os;
import sys;
from collections import defaultdict

import pysam;

#from nrlbio.filters_for_sam import *
#from nrlbio.chimera import Chimera, as_score
from nrlbio.samlib import ArWrapper
from nrlbio.generators import generator_doublesam

parser = argparse.ArgumentParser(description='produce chimeras from merged and already filtered sam file');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to merged sam/bam file, if the chimeric hits are present there consequently(one reference scenario). Path to two sam/bam files, if chimeric reads are splited between them(two reference scenario)");
parser.add_argument('--fixgap', nargs = '?', default = False, const=True, type = bool, help = "Can be set in two reference scenario, when already mapped left part was removed before the second mapping")
parser.add_argument('--oformat', nargs = '?', default = 'bed', choices=['bed', 'gff'], type = str, help = "format of the output file")
parser.add_argument('--collapsed', nargs = '?', default = False, const=True, type = bool, help = "This flag has to be set, if reads were collapsed");
args = parser.parse_args();

strand_conv = {True: '-', False: '+'}
if(args.oformat == 'bed'):
	bedformat = True;
else:
	bedformat = False;
	from nrlbio.pybedtools_extension import  construct_gff_interval

def arws2doublebed(arws, fixgap=False, bedformat=True):
	'''converts doublesam(alignment wrappers) into bed intervals
	
		arws list: list of alignment wrappers which represent mapping hits to fifferent loci for the sam reads
		fixgap bool: has to be used in two reference scenario, when already mapped left part was removed before the second mapping
	'''
	r = [];
	
	if(fixgap):
		gap = arws[1].qstart;
		#gap_seq = arws[1].aligned_read.query_sequence[:arws[1].qstart];
	else:	
		gap = arws[1].qstart - arws[0].qend;	
		
		#if(arws[0].aligned_read.is_reverse):
			#gap_seq = reverse_complement(arws[0].aligned_read.query_sequence)[arws[1].qstart:arws[0].qend]
		#else:
			#gap_seq = arws[0].aligned_read.query_sequence[arws[1].qstart:arws[0].qend]
		
	score = sum([x.score for x in arws])	
			
	if(bedformat):		
		for c, arw in enumerate(arws):
			r.append("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\n" % (arw.rname, arw.aligned_read.reference_start, arw.aligned_read.reference_end, "|".join((arw.qname, str(c))), arw.AS, strand_conv[arw.aligned_read.is_reverse], arw.qstart, arw.qend, score, gap));
	else:
		for c, arw in enumerate(arws): 
			r.append(str(construct_gff_interval(arw.rname, arw.aligned_read.reference_start, arw.aligned_read.reference_end,'ch', score=str(arw.AS), strand=strand_conv[arw.aligned_read.is_reverse], source='un', frame='.', attrs=[('ID', "|".join((arw.qname, str(c)))), ('gap', gap), ('chscore', score), ('qstart', arw.qstart), ('qend', arw.qend), ('n_uniq', arw.n_uniq) ])));
			
		
	return r;
		
		
		

if(len(args.path) == 1):
	samfile = pysam.Samfile(args.path[0])
	for segments in generator_doublesam(samfile):
		segchroms = [samfile.getrname(segment.tid) for segment in segments]
		arws = [ArWrapper(segment, segchrom, add_nr_tag=args.collapsed) for segment, segchrom in zip(segments, segchroms)]
		for interval_string in arws2doublebed(arws, args.fixgap, bedformat):
			sys.stdout.write(interval_string)
		
		
elif(len(args.path) == 2):
	chimeras = defaultdict(list);
	
	samfile = pysam.Samfile(args.path[0]);
	for segment in samfile.fetch(until_eof=True):
		segchrom = samfile.getrname(segment.tid) 
		chimeras[segment.query_name].append(ArWrapper(segment, segchrom, add_nr_tag=args.collapsed))
	samfile.close();
	
	samfile = pysam.Samfile(args.path[1]);
	for segment in samfile.fetch(until_eof=True):
		segchrom = samfile.getrname(segment.tid) 
		chimeras[segment.query_name].append(ArWrapper(segment, segchrom, add_nr_tag=args.collapsed))
	samfile.close();
	
	
	for arws in chimeras.values():
		if(len(arws)==2):
			for interval_string in arws2doublebed(arws, args.fixgap, bedformat):
				sys.stdout.write(interval_string)
else:
	sys.stderr.write('Incorrect number of input sam/bam files. It has to be 1 or 2\n')
