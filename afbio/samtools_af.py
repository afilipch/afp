# /usr/bin/python
'''collections of classes and functions to deal with sam/bam files'''

import sys;
from collections import namedtuple, Counter;
from math import log
from collections import defaultdict

from afbio import numerictools;
from afbio.sequencetools import entropy;


strand_conv = {True: '-', False: '+'}



def get_conversions(ar):
    '''function to get all mismatches/deletions/insertions for given pysam.AlignedRead

    ar pysam.AlignedRead: aligned read to be analyzed

    Return list of tuples. 1st element in each tuple is nucleotide/gap (string/None) in reference. 2st element in each tuple is nucleotide/gap (string/None) in query.
    '''	
    conversions = [];
    mmdict, deliter = intermediate_alignment(ar);
    ins_adjust = 0;	

    for i, j in fix_aligned_pairs(ar):
        if(i != None):
            if(j != None):
                if(i - ins_adjust in mmdict):
                    apair = (mmdict[i-ins_adjust], ar.query[i]); 
                else:
                    continue;
            else:
                ins_adjust += 1;
                apair = (None, ar.query[i])
        else:
            apair = (next(deliter), None)
        conversions.append(apair);	
            
    return conversions;


#__________________________________________________________________________________
#Local Auxillary functions
#__________________________________________________________________________________

def convert_positions_query(segment):
    '''For mappings on reverse strand, query start and end of the allignment are reversed and have to be reversed for chimeras. 
    This function returns \'true\' start and end position of an alignment on a query for both reverse and forward mappings.

    Returns int, int: \'true\' start and end position of an alignment on a query 
    '''
    if(segment.is_reverse):
        return segment.query_length-segment.query_alignment_end, segment.query_length-segment.query_alignment_start
    else:
        return segment.query_alignment_start, segment.query_alignment_end



#__________________________________________________________________________________
#Scoring functions
#__________________________________________________________________________________
def as_score(arw):
    return arw.AS;
    
def as_qstart_score(arw):
    qs = 2-log(arw.qstart+1);
    return arw.AS*(1+qs)

def as_qstart_entropy_score(arw):
    qs = 2-log(arw.qstart+1);
    e = (entropy(arw.aligned_read.query) - 1.5)
    if(e<0):
        e = e*5
    return arw.AS*(1+qs+e)

def as_qstart_rstart_score(arw):
    qs = 2-log(arw.qstart+1);
    rs = 1-log(arw.aligned_read.pos+1);
    return arw.AS*(1+qs+rs)

def as_qstart_rstart_entropy_score(arw):
    qs = 2-log(arw.qstart+1);
    rs = 1-log(arw.aligned_read.pos+1);
    e = (entropy(arw.aligned_read.query) - 1.5)
    if(e<0):
        e = e*5
    return arw.AS*(1+qs+rs+e)
#__________________________________________________________________________________	
	

#__________________________________________________________________________________
#Demultiplexing API
#__________________________________________________________________________________	

class ArWrapper(object):
    '''Wrapper for pysam.aligned read. Adds some additional fields
        
    Attributes:
        aligned_read pysam.aligned_read: read to wrap
        qname str: read id
        rname str: reference id
        AS float: alignment score
        control bool: if True, read comes from decoy
        conversions list of tuples: 1st element in each tuple is nucleotide/gap (string/None) in reference. 2st element in each tuple is nucleotide/gap (string/None) in query.
    '''	
        
    def __init__(self, aligned_read, rname, score_function=as_score, add_nr_tag = False, secondmate=False):
        self.aligned_read = aligned_read;
        self.qname = aligned_read.qname;
        self.rname = rname
        if(secondmate):
            self.strand = strand_conv[not aligned_read.is_reverse]
        else:
            self.strand = strand_conv[aligned_read.is_reverse]
        
        self.qstart, self.qend = convert_positions_query(aligned_read)
        
            
        self.AS = aligned_read.opt("AS")
        
        if(rname.split("_")[0] == "random"):
            self.control = True;
        else:
            self.control = False;
            
        if(add_nr_tag):
            self.n_uniq = int(self.qname.split("_")[-1][1:])
            self.aligned_read.tags = self.aligned_read.tags + [("NR", self.n_uniq)];
        else:
            self.n_uniq = 1; 
            
        self.score = score_function(self);
        
        
    def set_conv(self, from_, to):
        '''adds number of given type of conversions as a tag to the self.aligned_read
        #filtering API		
            from_ char: conversion from (from 'T' in PAR-CLIP)
            to char: conversion to (to 'C' in PAR-CLIP)
        '''
        self.conversions = get_conversions(self.aligned_read);
        
        conv_number = Counter(self.conversions)[(from_, to)];
        conv = "".join((from_, to));
        
        self.aligned_read.tags = self.aligned_read.tags + [(conv, conv_number)];
        
        
    def reassign_coordinates(self, reassign=True):
        '''reassigns coordinates coordinates to genomic ones. NOTE: If reassign=False, then only new attributes appear, without reassignment. This is necessary for downstream compatibility'''
        if(reassign):
            chrom, strand, start, stop = self.rname.split("|")[:4]
            self.chrom = chrom;
            self.strand = strand;
            start = int(start);
            stop = int(stop);
            
            if(strand == '+'):
                self.stop = start + int(self.aligned_read.aend)
                self.start = start + int(self.aligned_read.pos)
            else:
                self.start = stop - int(self.aligned_read.aend)
                self.stop = stop - int(self.aligned_read.pos)
            
        else:
            self.start = int(self.aligned_read.pos)
            self.stop = int(self.aligned_read.aend)
            self.chrom = self.rname
                
                
    def __str__(self):
        return "\t".join([str(x) for x in (self.qname, self.rname, self.qstart, self.qend, self.strand, self.score)]);



def demultiplex_read_hits(arwlist, bestdistance, backward=False):
    '''Demultiplex hits derived from the same read (choose the best ones on basis of its score)

        arwlist list: ArWrappers of the aligned reads(hits) derived from the same reads
        bestdistance int: minimal distance allowed between the best and the second best hit. If the actual distance is less, than hit will be assigned as nonunique
        
        Returns list: list of all valid (within bestdistance difference to the best hit) hits (ArWrapper objects)
    '''
    bestchoice = max(arwlist, key =  lambda x: x.score)
    return filter(lambda x: (bestchoice.score-x.score)<bestdistance, arwlist)



        
        
        
        
        
#__________________________________________________________________________________
#filtering API
#__________________________________________________________________________________
def get_attributes(ar, attributes):
    '''Converts aligned_read into list corresponding to the attributes provided. We need to do so, since some of attribute of the aligned_read are not accessible via getattr'''
    l = []
    for attr in attributes:
        if(hasattr(ar, attr)):
            l.append(getattr(ar, attr));
        else:
            l.append(ar.opt(attr));
    return l;
	
	
def get_attributes_masked(ar, attributes):
    '''Converts aligned_read into list corresponding to the attributes provided. We need to do so, since some of attribute of the aligned_read are not accessible via getattr'''
    l = []
    for attr in attributes:
        if(hasattr(ar, attr)):
            if(attr == 'qstart'):
                l.append(ar.qstart -  ar.seq.rfind('N'))
            else:	
                l.append(getattr(ar, attr));
        else:
            l.append(ar.opt(attr));
    return l;
	

	
def filter_generator(samfile, attributes, ga = get_attributes):
    '''Yields list of attributes corresponding to the aligned_reads in samfile provided. Each list will be used as entry in further filtering
        
        samfile pysam.Samfile: samfile to generate lists for further filtering
        attributes list: list of attributes important for filtering
    '''
    for aligned_read in samfile.fetch(until_eof=True):
        if(not aligned_read.is_unmapped):
            yield ga(aligned_read, attributes);
    samfile.close()		
            
        
def apply_filter(samfile, attributes, filter_, ga = get_attributes):
    '''Applies given filter to each entry(aligned) in the samfile
        
        samfile pysam.Samfile: samfile to generate lists for further filtering
        attributes list: list of attributes important for filtering. Important: it must be the same as attributes argument in filter_generator
        filter_ str: rule to filter list corresponding to each aligned_read
        
    Yields pysam.AlignedRead: sam entry passed the filtering	
    '''	
    for aligned_read in samfile.fetch(until_eof=True):
        if(not aligned_read.is_unmapped):
            x = ga(aligned_read, attributes)
            if(eval(filter_)):
                yield aligned_read
            else:
                pass;
    samfile.close();




#__________________________________________________________________________________	
#Miscellaneous
#__________________________________________________________________________________
def sort_segments_qstart(segments):
    segments.sort(key=lambda x: convert_positions_query(x)[0]);

            

