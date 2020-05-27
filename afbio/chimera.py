'''library contains classes and functions to deal with chimeric reads'''

import sys
from copy import copy
from itertools import combinations;

from Bio.Seq import reverse_complement


from afbio.numerictools import isinteger


class ChimeraException(Exception):
    pass

    
#__________________________________________________________________________________
#demultiplexing API

#__________________________________________________________________________________    
    
strand_conv = {True: '-', False: '+'}
    
def as_score(chimera):
    '''score function for Chimera based only on alignment score'''
    return sum(chimera.AS)




class Chimera(object):
	'''Chimera represents chimeric read based on two consecutive(with some gap or overlap) mappings for one read.
	
		aligned_reads list: element is pysam.AlignedRead. Two consecutive hits comprising the chimera. Must be in consecutive(from left to right) order
		rnames list: list of hits' reference names. The order must correspond to the order of aligned_reads list
		AS list: list of alignment scores corresponding to the aligned_reads
		self.gap int: gap(if negative, overlap) between hits(aligned_reads)
		self.score float: Chimera score. The bigger the score the more reliable is chimera
	'''	
		
		
	def __init__(self, arws):
		if(len(set([x.qname for x in arws])) != 1): 
			raise ChimeraException('Chimera cannot be made from aligned reads with different identifiers\nFollowing are given:\n%s\n' % "\n".join(["\t%s" % x for x.qname in ar.wrappers]))
			
		self.arws = arws;
		self.control = any([x.rname.split("_")[0] == "random" for x in self.arws])
		
		self.gap = arws[1].qstart - arws[0].qend;
		self.AS = sum([x.AS for x in arws]);
		self.coordinates = [];
		
		if(arws[0].aligned_read.is_reverse):
			self.gap_seq = reverse_complement(arws[0].aligned_read.query_sequence)[self.arws[1].qstart:self.arws[0].qend]
		else:
			self.gap_seq = arws[0].aligned_read.query_sequence[self.arws[1].qstart:self.arws[0].qend]
			
			
	def set_score(self, score):
		self.score = score
		
		
	def get_coordinates(self, reassign):
		
		if(reassign):
			for arw in self.arws:
				chrom, strand, start, end = arw.rname.split("|")[:4]
				if(strand == '+'):
					rend = int(start) + arw.aligned_read.reference_end
					rstart = int(start) + arw.aligned_read.reference_start
				else:
					rend = int(end) - arw.aligned_read.reference_start
					rstart = int(end) - arw.aligned_read.reference_end
				self.coordinates.append((chrom, strand, rstart, rend))
		else:
			for arw in self.arws:
				self.coordinates.append((arw.rname, strand_conv[arw.aligned_read.is_reverse], arw.aligned_read.reference_start, arw.aligned_read.reference_end))
		
		
	def __str__(self):
		'''converts Chimera in a bed-like entry(line). First six elements correspond to the first hit, Second six elements to the second one'''
		l = [];
		for arw in self.arws:
			l.append("%s\t%d\t%d\t%s\t%d\t%s\t" % (arw.rname, arw.aligned_read.pos, arw.aligned_read.aend, arw.qname, arw.AS, strand_conv[arw.aligned_read.is_reverse]))
			
		#gen_info = "%1.5f\t%d" % (self.score, self.gap);
		#l.append(gen_info)
	
		#for arw in self.arws:
			#l.append("\t%d\t%d\t%s" % (arw.qstart, arw.qend, arw.aligned_read.query[-1]))
		
		return "".join(l);	
		
		
	def __cmp__(self, other):
		return  cmp(self.score, other.score)
		
		
	def doublebed(self):
		'''Converts chimera in two bed entries with the same identifiers'''
		l = [];
		gen_info = "%1.5f\t%d" % (self.score, self.gap);
		
		if(self.gap_seq):
			gap_seq = self.gap_seq;
		else:
			gap_seq = str(self.gap)
		
		if(self.coordinates):
			for c, arw in enumerate(self.arws):
				l.append("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s" % (self.coordinates[c][0], self.coordinates[c][2], self.coordinates[c][3], "|".join((arw.qname, str(c))), arw.AS, self.coordinates[c][1], arw.qstart, arw.qend,  gen_info));			
		else:
			for c, arw in enumerate(self.arws):
				l.append("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s" % (arw.rname, arw.aligned_read.pos, arw.aligned_read.aend, "|".join((arw.qname, str(c))), arw.AS, strand_conv[arw.aligned_read.is_reverse], arw.qstart, arw.qend, gen_info));
			
		return "\n".join(l)
		
		
def arwlist2chimeras(arwlist, maxgap, maxoverlap):
	'''Compiles all possible chimeras from hits(pysam aligned reads) provided
	
		arlist iterable: element is pysam.AlignedRead. All hits of the initial read to compile into chimeras
		gap int: maximum gap allowed between two hits
		gap int: maximum gap allowed between two hits
		
		Returns list: all possible chimeras compiled from hits(pysam aligned reads)
	'''
	chimeras = []
	for a1, a2 in combinations(arwlist, 2):
		if(-maxoverlap <= a2.qstart - a1.qend <= maxgap):
			chimeras.append(Chimera([a1, a2]))
			continue
		if(-maxoverlap <= a1.qstart - a2.qend <= maxgap):
			chimeras.append(Chimera([a2, a1]))
	return chimeras;
	



	
def demultiplex(chimeras, bestdistance, backward=False):
	'''Demultiplex chimeric hits derived from the same read (choose the best ones on basis of its score)
	
		chimeras list: Chimeras of the aligned reads(hits) derived from the same reads
		bestdistance int: minimal distance allowed between the best and the second best chimera. If the actual distance is less, than  will chimera be assigned as nonunique
		
		Returns list: list of all valid (within bestdistance difference to the best chimera) chimera (Chimera objects)
	'''
	if(chimeras):
		bestchoice = max(chimeras, key =  lambda x: x.score)
		return filter(lambda x: (bestchoice.score-x.score)<bestdistance, chimeras)
	else:
		return []
			
		
		
		
#__________________________________________________________________________________
#filtering API
#__________________________________________________________________________________

def get_attributes(a, indices):
	'''Converts line in chimera file in a list ready to be passed to filtering'''
	#print indices
	#print a;
	for i in indices:
		if(isinstance(a[i], int)):
			pass;
		elif(isinstance(a[i], str) and isinteger(a[i])):
			a[i] = int(a[i]);
			#sys.stderr.write("%s\t%s\n" % (i, a[i]))
		else:
			a[i] = float(a[i]);
	return a;
	

	
def filter_generator(path, indices):
	'''Yields list of attributes corresponding to the chimera. Each list will be used as entry in further filtering
		
		path str:path to the file with chimeras
		indices list of integers: list of indices important for filtering
	'''
	with open(path) as f:
		for l in f:
			a = l.strip().split("\t");
			yield get_attributes(a, indices);
			
			
def filter_generator_doublebed(path, indices):
	'''USE FOR DOUBLEBED FORMAT. Yields list of attributes corresponding to the chimera. Each list will be used as entry in further filtering
		
		path str:path to the file with chimeras
		indices list of integers: list of indices important for filtering
	'''
	i = 1;
	with open(path) as f:
		for l in f:
			if(i % 2):
				a1 = l.strip().split("\t");
			else:
				a2 = l.strip().split("\t");
				a = a1[0:6] + a2[0:6] + a1[6:11] + a2[6:]
				yield get_attributes(a, indices);	
			i+=1	
			
			
def apply_filter(path, indices, filter_):
	'''Applies given filter to each entry(aligned) in chimeras
		
		path str:path to the file with chimeras
		indices list of integers: list of indices important for filtering
		filter_ str: rule to filter list corresponding to each chimera
		
	Yields str: ready to print representation of chimera
	'''	
	with open(path) as f:
		for l in f:
			a = l.strip().split("\t");
			x = get_attributes(a, indices)
			if(eval(filter_)):
				yield l;
			else:
				pass;
				

def apply_filter_doublebed(path, indices, filter_):
	'''USE FOR DOUBLEBED FORMAT. Applies given filter to each entry(aligned) in chimeras
		
		path str:path to the file with chimeras
		indices list of integers: list of indices important for filtering
		filter_ str: rule to filter list corresponding to each chimera
		
	Yields str: ready to print representation of chimera
	'''	
	for x in filter_generator_doublebed(path, indices):
			if(eval(filter_)):
				s1 = "\t".join([str(e) for e in x[0:6] + x[12:17]]);
				s2 = "\t".join([str(e) for e in x[6:12] + x[17:]]);
				yield "\n".join((s1, s2));
			else:
				pass;		



				
				
#testing section
if(__name__ == "__main__"):
	pass
