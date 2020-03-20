#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Selects the best hits from blastp output'''

import argparse
import os
import sys

import dominate
from dominate.tags import *
from dominate.util import raw
from Bio import SearchIO


parser = argparse.ArgumentParser(description='Selects the best hits from blastp output');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to blast output");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
#parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Path to the table with proteins from human genome atlas");
parser.add_argument('--identity', nargs = '?', default=0.5, type = float, help = "minimum protein identity");
args = parser.parse_args()



with open(args.css) as css:
    _style = css.read()
    
    
def sign_to_class(sign):
    if(sign == ' '):
        return 0;
    if(sign == '+'):
        return 1;
    return 2;
    
def convert_aln(hsp):
    res = []
    curstart = 0;
    alnann = hsp.aln_annotation['similarity'][:99]
    curclass = sign_to_class(alnann[0])
    
    for start, sign in enumerate(alnann[1:], start=1):
        tclass = sign_to_class(sign)
        if(tclass != curclass):
            res.append(( curstart, start, curclass))
            curclass = tclass;
            curstart = start;
    else:
        res.append(( curstart, len(alnann), curclass ))
        
    print(res)
    print(alnann[:99])
        

def create_hsp_html(hsp):
    #with open(os.path.join(args.outdir, "%s_%s.html" % gene.), 'w') as f:
    
    hsp.hit.id
    
    doc = dominate.document(title="%s VS %s" % (hsp.query.id, hsp.hit.id) )
    with doc.head:
        style(_style)
    with doc:
        p(strong("%s VS %s" % (hsp.query.id, hsp.hit.id)))
        br()
        with p():
            b("Query:")
            raw("\t%s\t%d-%d" % (hsp.query.id, hsp.query_start, hsp.query_end) )
        with p():
            b("Hit:")
            raw("\t%s\t%d-%d" % (hsp.hit.id, hsp.hit_start, hsp.hit_end) )
        #with p():
            #b("CDS:")
            #raw('<a href=%s target="_blank">\t%s(%s) %d-%d</a>' % (add_ucsc(gene, args.ucsc), gene.chrom, gene.strand, gene.start, gene.stop) )
            ##raw("\t%s(%s) %d-%d " % (gene.chrom, gene.strand, gene.start, gene.stop))
        #with p():
            #b("Alternative TSS:")
            #raw("\t%s" % gene.attrs["alt_tss"].replace(",", ", ") )
        #with p():
            #b("Annotation:")
            #raw("\t%s" % gene.attrs['annotation'] )
        #with p():
            #b("Function:")
            #raw("\t%s" % gene.attrs['function'] )
        #with p():
            #b("Phage:")
            #raw("\t%s" % PHAGE[gene.attrs.get('phage')] )
        #with p():
            #b("Protein seq:")
            #raw("\t%s" %  interval2seq(interval, genome))
            
    print(doc.render())
                
                
    


def find_identity(aln):
    return len([ 1 for x in zip(aln[0], aln[1]) if x[0] == x[1] ])/len(aln[0])
        

for qresult in SearchIO.parse(args.path, 'blast-xml'):
    for hits in qresult:
        for hsp in hits:
            identity = find_identity(hsp.aln)
            if(identity>args.identity):
                #create_hsp_html(hsp)
                convert_aln(hsp)
                #stop = 100
                #print(hsp.hit_range, hsp.hit.id, len(hsp.hit))
                #print(hsp.aln[0].seq[:stop]) 
                #print(hsp.aln_annotation['similarity'][:stop])
                #print(hsp.aln[1].seq[:stop])
                

            



                
