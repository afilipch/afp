#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Selects the best hits from blastp output'''

import argparse
import os
import sys

import numpy as np;
from pathlib import Path
from collections import defaultdict
import dominate
from dominate.tags import *
from dominate.util import raw
from Bio import SearchIO
import copy


parser = argparse.ArgumentParser(description='Selects the best hits from blastp output');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to blast output");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--identity', nargs = '?', default=0.5, type = float, help = "minimum protein identity");
parser.add_argument('--minscore', nargs = '?', default=40, type = float, help = "minimum required local (custom) alignment score");
args = parser.parse_args()


LOCAL_HTML_PATH = os.path.join(args.outdir, 'local')
LOCAL_GENES_PATH = os.path.join(args.outdir, 'genes')
Path(LOCAL_HTML_PATH).mkdir(parents = True, exist_ok=True)
Path(LOCAL_GENES_PATH).mkdir(parents = True, exist_ok=True)
#sys.exit()

colors = ['lightgray', 'indigo', 'forestgreen']
aln_signs = [' ', '+', '|']
aln_scores = [-1.5, 0.5, 1]



    
def sign_to_class(sign):
    if(sign == ' '):
        return 0;
    if(sign == '+'):
        return 1;
    return 2;


def extend_score(seed1, seed2):
    return min(seed1[2], seed2[2]) + aln_scores[0]*(seed2[0]-seed1[1])
    
def merge_seeds(seeds):
    if(len(seeds) == 1):
        return seeds[0]
    else:
        score = sum(x[2] for x in seeds) + sum([ aln_scores[0]*(x[1][0]-x[0][1]) for x in zip(seeds, seeds[1:]) ])
        return seeds[0][0], seeds[-1][1], score
    

def split_aln(hsp, minscore):
    alnann = hsp.aln_annotation['similarity']
    curstart = 0;
    curscore = 0;
    seeds = []
    for start, sign in enumerate(alnann):
        score = aln_scores[sign_to_class(sign)]
        if(score<0):
            if(curscore):
                seeds.append(( curstart, start, curscore))
                curscore = 0;
        else:
            if(not curscore):
                curstart = start
            curscore += score;
    else:
        if(curscore):
            seeds.append(( curstart, len(alnann), curscore ))
            
    
    #for s,e, score in seeds:
        #print(alnann[s:e], s, e, score)
        
    #print()
    newseeds = [];
    curseeds = [seeds[0]]
    for seed1, seed2 in zip(seeds, seeds[1:]):
        es = extend_score(seed1, seed2)
        if(es>0):
            curseeds.append(seed2)
        else:
            newseeds.append(merge_seeds(curseeds))
            curseeds = [seed2];
    else:                
        newseeds.append(merge_seeds(curseeds))
     
    
    return [x for x in newseeds if x[2]>minscore]
    


def convert_split(hsp, seed):
    alnann = hsp.aln_annotation['similarity'][seed[0]:seed[1]]
    identity = 1 - (alnann.count("+") + alnann.count(" "))/len(alnann)
    #print(alnann, seed, identity)
    res = []
    curstart = 0;
    curclass = sign_to_class(alnann[0])
    for start, sign in enumerate(alnann[1:], start=1):
        tclass = sign_to_class(sign)
        if(tclass != curclass):
            res.append(( curstart, start, curclass))
            curclass = tclass;
            curstart = start;
    else:
        res.append(( curstart, len(alnann), curclass ))
        
    return res, identity

 
def add_uniprot_ref(name): 
    return "https://www.uniprot.org/uniprot/%s" % name


def seed_to_html(hsp, seed, aln_marks, _style):
    h_start = hsp.hit_start + seed[0];
    h_end = hsp.hit_start + seed[1];
    q_start = hsp.query_start + seed[0];
    q_end = hsp.query_start + seed[1];
    _ , uniprot, name = hsp.hit.id.split("|")
    
    doc = dominate.document(title="Query: %s Hit: %s Hit Range: %d-%d" % (hsp.query.id, name, h_start, h_end) )
    with doc.head:
        style(_style)   
    with doc:
        with p():
            b("Query: ")
            raw(hsp.query.id)
            b("Range: ")
            raw("%d-%d" % (q_start, q_end) )           
        with p():
            b("Hit: ")
            raw('<a href=%s target="_blank">%s</a>' % (add_uniprot_ref(uniprot), name) )
            b("Range: ")
            raw("%d-%d" % (h_start, h_end) ) 
        with p():
            b("Score: ")
            raw("%1.1f" % seed[2])
        with p():
            b("Identity: ")
            raw("\t%1.2f%%" % (100*identity))
        with p():
            b("Alignment:")
            seq1 = str(hsp.aln[0].seq).upper()[seed[0]:seed[1]]
            seq2 = str(hsp.aln[1].seq).upper()[seed[0]:seed[1]]
            with p():
                ugly_fonts = [[], []]
                aln_line = ''
                for start, stop, tclass in aln_marks:
                    s = '<font color=%s>%s</font>' % (colors[tclass], seq1[start:stop])
                    ugly_fonts[0].append(s);
                    s = '<font color=%s>%s</font>' % (colors[tclass], seq2[start:stop])
                    ugly_fonts[1].append(s);
                    aln_line += aln_signs[tclass]*(stop-start)
                    
                raw(''.join(ugly_fonts[0]))
                br()
                font(aln_line, color='black')
                br()
                raw(''.join(ugly_fonts[1]))
                
    fname = "%s_%s_%d_%d.html" % (hsp.query.id, name, h_start, h_end)
    fpath = os.path.join(LOCAL_HTML_PATH, fname)
    with open(fpath, 'w') as f:        
        f.write(doc.render())
    return (hsp.query.id, hsp.hit.id, fname, q_start, q_end, h_start, h_end, identity, seed[2])
        
        
########### EXECUTION SECTION
      
with open(args.css) as css:
    _style = css.read()
    
res_list = []
try:
    for qresult in SearchIO.parse(args.path, 'blast-xml'):
        #print('bu')
        #print("bu")
        for hits in qresult:
            for hsp_l in hits:
                for hsp in hsp_l:
                    seeds = split_aln(hsp, args.minscore)
                    for seed in seeds:
                        aln_marks, identity = convert_split(hsp, seed)
                        res = seed_to_html(hsp, seed, aln_marks, _style)
                        res_list.append(res)
except:
    sys.stderr.write("\nerror of xml parser skipped\n%d items were already collected\n" % len(res_list))
  
#print(len(res_list))

res_dict = defaultdict(list)
for el in res_list:
    res_dict[el[1]].append(el)
    
gene_level_res = []
header = "Query", "Hit", "Link", "q_start", "q_end", "h_start", "h_end", "Identity [%]", "Score"
for hit, local_list in res_dict.items():
    name_ = hit.split("|")[-1]
    gene_level_res.append(( hit, len(local_list), sum([x[-1] for x in local_list]) ))
    local_list.sort(key = lambda x: x[-1], reverse = True)
    with open(os.path.join(LOCAL_GENES_PATH, '%s.html' % name_), 'w') as f:
        doc = dominate.document(title="%s: similarity blast ppredictions" % os.path.basename(args.outdir))
        with doc.head:
            style(_style)
        with doc:
            with table(id = "myTable") as _table:
                _tr = tr()
                _tr.add([td(x) for x in header])
                for query, hit, fname, q_start, q_end, h_start, h_end, identity, score in local_list:
                    _ , uniprot, name = hit.split("|")
                    with tr():
                        td(query)
                        td(raw('<a href=%s target="_blank">%s</a>' % (add_uniprot_ref(uniprot), name) ))
                        td(raw('<a href=%s target="_blank">link</a>' % os.path.join('../local', fname) ))
                        td(q_start)
                        td(q_end)
                        td(h_start)
                        td(h_end)
                        td("%1.2f" % (100*identity))
                        td("%1.1f" % score)
        f.write(doc.render())
        
        

header = "gene name", "link to hits", "number of hits", "total score"
gene_level_res.sort(key = lambda x: x[-1], reverse = True)      
doc = dominate.document(title="%s: Human genes with counterparts" % os.path.basename(args.outdir))
with doc.head:
    style(_style)
with doc:
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([td(x) for x in header])
        for hit, num, score in gene_level_res:
            _ , uniprot, name = hit.split("|")
            with tr():
                td(raw('<a href=%s target="_blank">%s</a>' % (add_uniprot_ref(uniprot), name) ))
                td(raw('<a href=%s target="_blank">link</a>' % os.path.join('genes', "%s.html" % name) ))
                td(num)
                td("%1.1f" % score)
                
with open(os.path.join(args.outdir, 'blastp_genes.html'), 'w') as f:
    f.write(doc.render())
 
                

 
 
 
res_list.sort(key = lambda x: x[-1], reverse = True)
with open(os.path.join(args.outdir, 'blastp_local_hits.tsv'), 'w') as f:
    for res in res_list:
        f.write("\t".join(list(res[:2]) + [str(x) for x in res[3:]]) + "\n")
                

            


### OLD functions

#def find_identity(aln):
    #return len([ 1 for x in zip(aln[0], aln[1]) if x[0] == x[1] ])/len(aln[0])

#def convert_aln(hsp, size = 100):
    #total_ann = hsp.aln_annotation['similarity']
    #total_res = []
    #for div_start in range(0, len(total_ann), size):
        #alnann = total_ann[div_start:div_start+size]
    
        #res = []
        #curstart = 0;
        #curclass = sign_to_class(alnann[0])
        
        #for start, sign in enumerate(alnann[1:], start=1):
            #tclass = sign_to_class(sign)
            #if(tclass != curclass):
                #res.append(( curstart, start, curclass))
                #curclass = tclass;
                #curstart = start;
        #else:
            #res.append(( curstart, len(alnann), curclass ))
            
        #total_res.append(( div_start, copy.copy(res) ))
        ##print(res)
        ##print()
    #return total_res;
    
    
    
#def create_hsp_html(hsp, aln_marks, identity, fpath, evalue):

    #doc = dominate.document(title="%s VS %s" % (hsp.query.id, hsp.hit.id) )
    #with doc.head:
        #style(_style)
    #with doc:
        #p(strong("%s VS %s" % (hsp.query.id, hsp.hit.id)))
        #br()
        #with p():
            #b("Query:")
            #raw("\t%s\t%d-%d" % (hsp.query.id, hsp.query_start, hsp.query_end) )
        #with p():
            #b("Hit:")
            #raw('<a href=%s target="_blank">\t%s\t%d-%d</a>' % (add_uniprot_ref(hsp.hit.id), hsp.hit.id, hsp.hit_start, hsp.hit_end) )
        #with p():
            #b("Identity:")
            #raw("\t%1.2f%%" % (100*identity))
        #with p():
            #b("e-value:")
            #raw(str(evalue))
            
        #with p():
            #b("Alignment:")
            #seq1 = str(hsp.aln[0].seq).upper()
            #seq2 = str(hsp.aln[1].seq).upper()
            #for div_start, l1 in aln_marks:
                #with p():
                    #b("start: %d\t" % div_start)
                    #br()
                    #ugly_fonts = [[], []]
                    #aln_line = ''
                    #for start, stop, tclass in l1:
                        #s = '<font color=%s>%s</font>' % (colors[tclass], seq1[div_start+start:div_start + stop])
                        #ugly_fonts[0].append(s);
                        #s = '<font color=%s>%s</font>' % (colors[tclass], seq2[div_start+start:div_start + stop])
                        #ugly_fonts[1].append(s);
                        #aln_line += aln_signs[tclass]*(stop-start)
                        
                    #raw(''.join(ugly_fonts[0]))
                    #br()
                    #font(aln_line, color='black')
                    #br()
                    #raw(''.join(ugly_fonts[1]))
                        ##font(seq[div_start+start:div_start + stop], color = colors[tclass]);
    #with open(fpath, 'w') as f:        
        #f.write(doc.render())
    
    
    
    
#res_list = []
#for qresult in SearchIO.parse(args.path, 'blast-xml'):
    #for hits in qresult:
        #for hsp_l in hits:
            #for hsp in hsp_l:
                #identity = find_identity(hsp.aln)
                #if(identity>args.identity):
                    #aln_marks = convert_aln(hsp)
                    #fname = "%s_%s_%d_%d.html" % (hsp.query.id, hsp.hit.id, hsp.query_start, hsp.query_end)
                    #fpath = os.path.join(LOCAL_HTML_PATH, fname)
                    #create_hsp_html(hsp, aln_marks, identity, fpath, hsp_l.evalue)
                    #res_list.append(( hsp, identity, fname, hsp_l.evalue )) 
 


#header = "Query", "Hit", "Link", "q_start", "q_end", "h_start", "h_end", "identity [%]", "e-value"
#_title = "%s: similarity blast ppredictions" % os.path.basename(args.outdir)
#doc = dominate.document(title=_title)

#with open(args.css) as f:
    #_style = f.read()

#with doc.head:
    #style(_style)

#with doc:
    #p(strong(_title))
    #with table(id = "myTable") as _table:
        #_tr = tr()
        #_tr.add([td(x) for x in header])
        #for hsp, identity, fname, evalue in res_list:
            #with tr():
                #td(hsp.query.id)
                #td(raw('<a href=%s target="_blank">%s</a>' % (add_uniprot_ref(hsp.hit.id), hsp.hit.id) ))
                #td(raw('<a href=%s target="_blank">link</a>' % os.path.join('local', fname) ))
                #td(hsp.query_start)
                #td(hsp.query_end)
                #td(hsp.hit_start)
                #td(hsp.hit_end)
                #td("%1.2f" % (100*identity))
                #td(evalue)
 
                
#with open(os.path.join(args.outdir, 'report.html'), 'w') as f:
    #f.write(doc.render())
                
