# /usr/bin/python
'''generators for different kinds of input files. They may be used to allow multiprocessing and save memory footprint.
usually you should provide path to input file and size of elements(not lines) in buffer'''
import sys;
import re;
import os;
import itertools;
from afbio import sequencetools
from collections import *


fastq = namedtuple('fastq', 'id, seq, sign, qual')

def generator_fastq(path, take = ['id', 'seq', 'sign', 'qual'], reverse = False, shuffle = False):
    '''yields fastq object(collections.namedtuple based) from fastq file provided

    path str: path to fastq file
    take list: list of fastq entry attributes to pass to yielded fastq object. Allows to save some memory. Entries have to be: "id", "seq", "sign", "qual"
    reverse bool: if True fastq sequence and quality string will be reversed
    shuffle bool: if True fastq sequence will be reversed

    Yields fastq: fastq object(collections.namedtuple based) 
    '''
    my_take = {"id": 0, "seq": 1, "sign": 2, "qual": 3}
    line_set = [];
    
    with open(path) as f:
        for line in f:  
            line_set.append(line.strip());
            if(len(line_set) == 4):
                appendix = [None]*4;
                for el in take:
                    if(reverse and el in ["seq", "qual"]):
                        appendix[my_take[el]] = line_set[my_take[el]][::-1];
                    elif(shuffle and el == "seq"):
                        appendix[my_take[el]] = sequencetools.shuffle_string(line_set[my_take[el]])
                    else:
                        appendix[my_take[el]] = line_set[my_take[el]]
                yield fastq(*appendix)
                line_set = [];

	
def generator_maf(path, aligned_species=None):
    '''reads maf stitched file, taking into account interactions id. Yields conservation.Maf object corresponding to each entry in maf file
    
    path string: path to maf stitched file
    aligned_species list: names of aligned species(genome suystems: ce6, mm9, etc.) to be considered in further analysis. If not given takes all species
    
    Yields conservation.Maf: object corresponding to each entry in maf file
    '''	
    
    from conservation import Maf;
    
    nonaligned = re.compile('-+$')
    alignment = OrderedDict()
    with open(path) as f:
        for l in f:
            l = l.strip();  
            if(l):
                if(l[0] == ">"):
                    if(len(l.split(":")) == 3):
                        header = l
                        mflag = True;
                    else:
                        specie = l[1:];
                        mflag = False;
                elif(aligned_species == None or specie in aligned_species):
                    if(mflag):
                        refseq = l.upper();
                    elif(not nonaligned.match(l)):
                        alignment[specie] = l.upper();
            else:
                yield Maf(header, alignment, refseq)
        try:
            yield Maf(header, alignment)
        except:
            pass
            
		
def generator_doublebed(path):
    '''Reads \'double\' bed file. Yields consecutive pairs of bed/gff intervals
    
    path string: path to doublebed file
    
    Yields list: 2-element list of BedTool intervals, representing chimera or interaction.
    '''

    from pybedtools import BedTool;	
    
    doublebed  = [];
    for c, interval in enumerate(BedTool(path)):
        doublebed.append(interval);
        if(c % 2):
            yield doublebed;
            doublebed = [];
                    

def generator_doublesam(samfile):
    '''Reads \'double\' sam file(chimeras). Yields consecutive pairs of sam segments
    
    samfile pysam.Samfile: sam file python wrapper
    
    Yields list: 2-element list of pysam segments, representing chimera
    '''
    
    doublesam  = [];
    for c, segment in enumerate(samfile.fetch(until_eof=True)):
        doublesam.append(segment);
        if(c % 2):
            seq1, seq2 = doublesam
            if(seq1.query_name != seq2.query_name):
                sys.stderr.write("WARNING: Consequtive paired reads have different names:\t%s\t%s\n" % (seq1.query_name, seq2.query_name))
            yield doublesam;
            doublesam = [];
            

                    

def generator_seqrecord(paths, ftype):
    '''Generates Bio.SeqRecord.SeqRecord objects from multiple files

    paths list: list of paths to files with seqrecords(genbank, fastq, fasta, ets.)
    ftype str: format of files with seqrecords('genbank', 'fastq', 'fasta', ets.)
        
    Yields Bio.SeqRecord.SeqRecord
    '''
    from Bio import SeqIO;	
    for path in paths:
        for seqrecord in SeqIO.parse(path, ftype):
            yield seqrecord;



def generator_seq(paths, ftype):
    '''Generates sequences from fasta/fastq/genbank files. NOTE: sequences will be converted to upper case

        paths list: list of paths to files with seqrecords(genbank, fastq, fasta, ets.)
        ftype str: format of files with seqrecords('genbank', 'fastq', 'fasta', ets.)
            
    Yiels str: nucleotide sequence
    '''
    from Bio import SeqIO;	
    for path in paths:
        for seqrecord in SeqIO.parse(path, ftype):
            yield str(seqrecord.seq.upper());
                
			
			
def generator_mirna(paths, seed_start=1, seed_stop=7):
    '''Generates miRNAs from given fasta files
    
        paths list: list of paths to fasta files with miRNAs
        seed_start int: start postion of mirna seed, 0-based and inclusive
        seed_stop int: stop postion of mirna seed, 0-based and exclusive
            
    Yields nrlbio.Mirna object
    '''
    from Bio import SeqIO;
    from nrlbio.mirna import Mirna
    for path in paths:
        for seqrecord in SeqIO.parse(path, 'fasta'):
            yield Mirna(seqrecord.id, str(seqrecord.seq.upper()), seed_start = seed_start, seed_stop = seed_stop)


def generator_segments(path, key_score= lambda x: x.AS, add_nr_tag=False, secondmate=False, converted=False):
    '''Generates all mapping hits for a single read from a given sam/bam file (path to it)'''

    import pysam;
    if(converted):
        from nrlbio.samlib import BackwardWrapper as CWrapper
    else:
        from nrlbio.samlib import ArWrapper as CWrapper
    
    samfile = pysam.Samfile(path)
    arwlist = [];
    current_name = '';
    
    for segment in samfile.fetch(until_eof=True):
        if(not segment.is_unmapped):
            name = samfile.getrname(segment.tid)
            arw = CWrapper(segment, rname, score_function=key_score, add_nr_tag=add_nr_tag, secondmate=secondmate)
            
            if(current_name and current_name != arw.qname):
                yield arwlist
                arwlist = [arw];

            else:
                arwlist.append(arw);
        else:
            if(current_name):
                yield arwlist;
                arwlist = [];
                    
        current_name = arw.qname;
                            
    else:
        yield arwlist
            


def targets_generator(consfasta):
    '''generates aligned blocks from fasta file(got from extract_maf.py) as python dictionaries'''
    from Bio import SeqIO
    
    sequences = {}
    for seqrecord in SeqIO.parse(consfasta, 'fasta'):
        a = seqrecord.id.split("|");
        if(len(a) == 6):
            if(sequences):
                yield name, mirid, sequences
            sequences = {}
            name, mirid = a[4], a[5]
        else:
            sequences[seqrecord.id] = str(seqrecord.seq.upper()).replace('U', 'T')
    else:
        yield name, mirid, sequences




def grouper(iterable, n):
    it = iter(iterable);
    arr = [];
    while(it):
        for i in range(n):
            try:
                arr.append(next(it));
            except:
                yield arr;
                return
        yield arr;
        arr = []
        
def get_only_files(folder):
    return [os.path.join(folder, f) for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) ]


def generator_paired_mappings(samfile):
    curname = '';
    readmappings = [];
    for seq1, seq2 in generator_doublesam(samfile):
        if(seq1.query_name == curname):
            readmappings.append((seq1, seq2))
        else:
            yield ([x for x in readmappings if x[0].is_proper_pair])
            readmappings = [(seq1, seq2)]
            curname = seq1.query_name
    else:
        yield ([x for x in readmappings if x[0].is_proper_pair])
        
        
        
def generator_single_mappings(samfile):
    curname = '';
    readmappings = [];
    for seq in samfile.fetch(until_eof=True):
        if(seq.query_name == curname):
            readmappings.append(seq)
        else:
            yield readmappings
            readmappings = [seq]
            curname = seq.query_name
    else:
        yield readmappings
    
    
       
    
    
    
                            
#testing
if(__name__ == "__main__"):
    pass;









































