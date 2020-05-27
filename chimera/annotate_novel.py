#! /usr/bin/python
'''Script tries to de novo annotate chimeras. To achieve it, linear and circular splice junctions are identified by applying following criteria:
1) No more than 100 kilobase distance between chimeric hits
2) GT/AG signal flanking the splice sites
3) Unambiguous chimeric breakpoint detection
''' 
import argparse
import sys;
from collections import defaultdict, Counter

from Bio import SeqIO


from afbio.generators import generator_doublebed
from afbio.numerictools import distance as find_distance
from afbio.sequencetools import seqrecord2seq
from afbio.pybedtools_af import construct_gff_interval, bed2gff




parser = argparse.ArgumentParser(description='Script tries to de novo annotate chimeras. To achieve it, linear and circular splice junctions are identified by applying following criteria:\n1) No more than 100 kilobase distance between chimeric hits\n2) GT/AG signal flanking the splice sites\n3) Unambiguous chimeric breakpoint detection', formatter_class = argparse.RawTextHelpFormatter);

parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the chimeras, double bed/gff file");
parser.add_argument('--distance', nargs = '?', default = 100000, type = int, help = "max allowed disctance between chimeric hits");
parser.add_argument('--reference', nargs = '?', required = True, type = str, help = "path to the reference(genome) to extract sequences from");
parser.add_argument('--stranded', nargs = '?', default = False, const=True, type = bool, help = "Should be set if the sequencing data are stranded")
parser.add_argument('--reverse', nargs = '?', default = False, const=True, type = bool, help = "Should be set if the sequencing data are reversed (For example - second mate reads)")
args = parser.parse_args();


def check_distance(i1, i2, distance):
    return i1.chrom==i2.chrom and i1.strand==i2.strand and find_distance((i1.start, i1.end), (i2.start, i2.end))<distance;



def check_splice_site(i1, i2, seqrecord, stranded, reverse):
    gap = int(i1.attrs['gap'])
    if(i1.strand == '+'):
        flank5 = seqrecord2seq(seqrecord, i1.end+gap, i1.end+2, strand='+')
        flank3 = seqrecord2seq(seqrecord, i2.start-2, i2.start-gap, strand='+') 
        #if len(flank3)==2:
            #print flank3;
    elif(i1.strand == '-'):
        flank5 = seqrecord2seq(seqrecord, i1.start-2, i1.start-gap, strand='-') 
        flank3 = seqrecord2seq(seqrecord, i2.end+gap, i2.end+2, strand='-')
        #if len(flank3)==2:
            #print flank3;
    else:
        sys.stderr.write("Strand is not defined assumed to be plus\n");
        
    breakpoints = [];
    antisense = False;
    for b in range(abs(gap)+1):
        if(reverse):
            if(flank5[b:b+2] == "CT" and flank3[b:b+2] == "AC"):
                breakpoints.append(b);
                antisense = True;
        else:		
            if(flank5[b:b+2] == "GT" and flank3[b:b+2] == "AG"):
                breakpoints.append(b);
            if(not stranded):
                if(flank5[b:b+2] == "CT" and flank3[b:b+2] == "AC"):
                    breakpoints.append(b);
                    antisense = True;
    return breakpoints, antisense;





def clarify(i1, i2, breakpoint, antisense):
    gap = int(i1.attrs['gap']);
    #clarify breakpoint
    #clarify strandness
    if(antisense):
        if(i1.strand=='+'):
            cs1 = construct_gff_interval(chrom=i2.chrom, start=i2.start+breakpoint, stop=i2.stop, feature='ch', score=i2.score, strand='-', source='.', frame='.', attrs=i2.attrs.items())
            cs2 = construct_gff_interval(chrom=i1.chrom, start=i1.start, stop=i1.stop+gap+breakpoint, feature='ch', score=i1.score, strand='-', source='.', frame='.', attrs=i1.attrs.items())
        else:
            cs1 = construct_gff_interval(chrom=i2.chrom, start=i2.start, stop=i2.stop-breakpoint, feature='ch', score=i2.score, strand='+', source='.', frame='.', attrs=i2.attrs.items())
            cs2 = construct_gff_interval(chrom=i1.chrom, start=i1.start-gap-breakpoint, stop=i1.stop, feature='ch', score=i1.score, strand='+', source='.', frame='.', attrs=i1.attrs.items())
    else:
        if(i1.strand=='+'):
            cs1 = construct_gff_interval(chrom=i1.chrom, start=i1.start, stop=i1.stop+gap+breakpoint, feature='ch', score=i1.score, strand=i1.strand, source='.', frame='.', attrs=i1.attrs.items())
            cs2 = construct_gff_interval(chrom=i2.chrom, start=i2.start+breakpoint, stop=i2.stop, feature='ch', score=i2.score, strand=i2.strand, source='.', frame='.', attrs=i2.attrs.items())
        else:	
            cs1 = construct_gff_interval(chrom=i1.chrom, start=i1.start-gap-breakpoint, stop=i1.stop, feature='ch', score=i1.score, strand=i1.strand, source='.', frame='.', attrs=i1.attrs.items())
            cs2 = construct_gff_interval(chrom=i2.chrom, start=i2.start, stop=i2.stop-breakpoint, feature='ch', score=i2.score, strand=i2.strand, source='.', frame='.', attrs=i2.attrs.items())
        
    return cs1, cs2;


def assign_splice_type(i1, i2):
    if(i1.strand=='+'):
        if(i1.start<i2.start):
            return 'lsj';
        else:
            return 'csj';
    else:
        if(i1.start>i2.start):
            return 'lsj';
        else:
            return 'csj';

		
	
	
	
reference = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

total = 0;
passed = 0;
uniqbreakpoints = [];
stypes = defaultdict(int);
for i1, i2 in generator_doublebed(args.path):
    #sys.stderr.write(str(i1))
    if(i1.file_type=='bed'):
        attrs1 = (('qstart', i1[6]), ('qend', i1[7]), ('chscore', i1[8]), ('gap', i1[9]))
        gi1 = bed2gff(i1, feature='ch', attrs=attrs1)
        attrs2 = (('qstart', i2[6]), ('qend', i2[7]), ('chscore', i2[8]), ('gap', i2[9]))
        gi2 = bed2gff(i2, feature='ch', attrs=attrs2);
    else:
        gi1 = i1;
        gi2 = i2;	

    total += 1;
    if(check_distance(gi1, gi2, args.distance)):
        seqrecord = reference[gi1.chrom];
        breakpoints, antisense = check_splice_site(gi1, gi2, seqrecord, args.stranded, args.reverse);
        uniqbreakpoints.append(len(breakpoints))
        if(len(breakpoints)==1):
            cs1, cs2 = clarify(gi1, gi2, breakpoints[0], antisense)
            splice_type = assign_splice_type(cs1, cs2)
            passed += 1;
        else:
            splice_type = 'intra'
            cs1, cs2 = gi1, gi2
    else:
        splice_type = 'inter'
        cs1, cs2 = gi1, gi2
        
    cs1.attrs['ntype'] = splice_type
    cs2.attrs['ntype'] = splice_type
    sys.stdout.write("%s%s" % (cs1, cs2));	
    stypes[splice_type]+=1;
    
    


num_splice_sites = Counter(uniqbreakpoints);
sys.stderr.write("total chimeras: %d\npassed chimeras: %d\nfraction passed %1.5f\n\n" % (total, passed, float(passed)/total));
sys.stderr.write("Ambiguity	in breackpoint detection\nnum canonical splice sites\tnum chimeras\n");
for k in sorted(num_splice_sites.keys()):
    sys.stderr.write("%d\t%d\n" % (k, num_splice_sites[k]));
sys.stderr.write("\nSplice site type\tnumber of reads\n");
for k, v in stypes.items():
    sys.stderr.write("%s\t%d\n" % (k, v));






