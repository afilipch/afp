#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Selects the best hits from blastp output'''

import argparse
import os
import sys

from pathlib import Path
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
args = parser.parse_args()


local_html_path = os.path.join(args.outdir, 'local')
Path(local_html_path).mkdir(parents = True, exist_ok=True)

#sys.exit()

colors = ['lightgray', 'indigo', 'forestgreen']
aln_signs = [' ', '+', '|']

with open(args.css) as css:
    _style = css.read()

def add_uniprot_ref(name): 
    return "https://www.uniprot.org/uniprot/%s" % name.split("|")[1]
    
def sign_to_class(sign):
    if(sign == ' '):
        return 0;
    if(sign == '+'):
        return 1;
    return 2;
    
def convert_aln(hsp, size = 100):
    total_ann = hsp.aln_annotation['similarity']
    total_res = []
    for div_start in range(0, len(total_ann), size):
        alnann = total_ann[div_start:div_start+size]
    
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
            
        total_res.append(( div_start, copy.copy(res) ))
        #print(res)
        #print()
    return total_res;
        
    
        

def create_hsp_html(hsp, aln_marks, identity, fpath, evalue):

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
            raw('<a href=%s target="_blank">\t%s\t%d-%d</a>' % (add_uniprot_ref(hsp.hit.id), hsp.hit.id, hsp.hit_start, hsp.hit_end) )
        with p():
            b("Identity:")
            raw("\t%1.2f%%" % (100*identity))
        with p():
            b("e-value:")
            raw(str(evalue))
            
        with p():
            b("Alignment:")
            seq1 = str(hsp.aln[0].seq).upper()
            seq2 = str(hsp.aln[1].seq).upper()
            for div_start, l1 in aln_marks:
                with p():
                    b("start: %d\t" % div_start)
                    br()
                    ugly_fonts = [[], []]
                    aln_line = ''
                    for start, stop, tclass in l1:
                        s = '<font color=%s>%s</font>' % (colors[tclass], seq1[div_start+start:div_start + stop])
                        ugly_fonts[0].append(s);
                        s = '<font color=%s>%s</font>' % (colors[tclass], seq2[div_start+start:div_start + stop])
                        ugly_fonts[1].append(s);
                        aln_line += aln_signs[tclass]*(stop-start)
                        
                    raw(''.join(ugly_fonts[0]))
                    br()
                    font(aln_line, color='black')
                    br()
                    raw(''.join(ugly_fonts[1]))
                        #font(seq[div_start+start:div_start + stop], color = colors[tclass]);
    with open(fpath, 'w') as f:        
        f.write(doc.render())
                
                
    


def find_identity(aln):
    return len([ 1 for x in zip(aln[0], aln[1]) if x[0] == x[1] ])/len(aln[0])
      
      
res_list = []
for qresult in SearchIO.parse(args.path, 'blast-xml'):
    for hits in qresult:
        for hsp_l in hits:
            for hsp in hsp_l:
                identity = find_identity(hsp.aln)
                if(identity>args.identity):
                    aln_marks = convert_aln(hsp)
                    fname = "%s_%s_%d_%d.html" % (hsp.query.id, hsp.hit.id, hsp.query_start, hsp.query_end)
                    fpath = os.path.join(local_html_path, fname)
                    create_hsp_html(hsp, aln_marks, identity, fpath, hsp_l.evalue)
                    res_list.append(( hsp, identity, fname, hsp_l.evalue )) 
 


header = "Query", "Hit", "Link", "q_start", "q_end", "h_start", "h_end", "identity [%]", "e-value"
_title = "%s: similarity blast ppredictions" % os.path.basename(args.outdir)
doc = dominate.document(title=_title)

with open(args.css) as f:
    _style = f.read()

with doc.head:
    style(_style)

with doc:
    p(strong(_title))
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([td(x) for x in header])
        for hsp, identity, fname, evalue in res_list:
            with tr():
                td(hsp.query.id)
                td(raw('<a href=%s target="_blank">%s</a>' % (add_uniprot_ref(hsp.hit.id), hsp.hit.id) ))
                td(raw('<a href=%s target="_blank">link</a>' % os.path.join('local', fname) ))
                td(hsp.query_start)
                td(hsp.query_end)
                td(hsp.hit_start)
                td(hsp.hit_end)
                td("%1.2f" % (100*identity))
                td(evalue)
 
                
with open(os.path.join(args.outdir, 'report.html'), 'w') as f:
    f.write(doc.render())
                

            



                
