import sys
import os
import dominate
import argparse
from dominate.tags import *
from dominate.util import raw
from pybedtools import BedTool

from afbio.html.methods import add_ucsc

parser = argparse.ArgumentParser(description='Creates html table for the annotated ChIP binding peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the assigned differential table, tsv format");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the gene annotation");
parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
parser.add_argument('--name', nargs = '?', default="unknown", type = str, help = "Name of the sample");
parser.add_argument('--top', nargs = '?', default=300, type = int, help = "Shows only [top] peaks");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory (tsv and html files will be written there)")
parser.add_argument('--genes_dir', nargs = '?', required=True, type = str, help = "Path to the directory with genes html reports")
args = parser.parse_args();





    

_title = "Differential gene expression of: %s" % args.name

data = []
with open(args.path) as f:
    header = next(f).strip().split("\t")
    for l in f:
        data.append(l.strip().split("\t"))
        
data.sort(key = lambda x: float(x[2]), reverse=True)
gene2annotation = dict([ (x.attrs['ID'], x) for x in BedTool(args.annotation)])


def annotate(d, gene2annotation):
    gene = gene2annotation[d[0]]
    return [gene.attrs['genesymbol']] + d[1:], gene
    






doc = dominate.document(title=_title)

with open(args.css) as f:
    _style = f.read()
    
with open(args.js) as f:
    plain_script = f.read()



with doc.head:
    style(_style)


header.insert(1, 'UCSC')
with doc:
    p(strong(_title))
    input(type="text", id="inp1", onkeyup='my_search(1, "inp1", "myTable")', placeholder="Search for a genesymbol..")
    #input(type="text", id="inp2", onkeyup='my_search(2, "inp2", "myTable")', placeholder="Search for a geneID..")
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x) for x in header])
        for d in data:
            newd, gene = annotate(d, gene2annotation)
            with tr():
                td(raw('<a href=%s target="_blank">%s</a>' % (os.path.join(args.genes_dir, gene.name + ".html"), newd[0]) ))
                td(raw('<a href=%s target="_blank">link</a>' % add_ucsc(gene, args.ucsc) ))
                for el in newd[1:]:
                    a = el.split(",")
                    if(len(a)==1):
                        td(el)
                    else:
                        with td():
                            for aa in a:
                                p(aa)
                
    _script = script(raw(plain_script), type='text/javascript')
    
        

with open(os.path.join(args.outdir, args.name + ".html" ), 'w') as f:
    #print(doc.render())
    f.write(doc.render());
    

with open(os.path.join(args.outdir, args.name + ".tsv"), 'w') as f:
    f.write("%s\n" % "\t".join(header[1:]));
    for d in data:
        newd, gene = annotate(d, gene2annotation)
        f.write("%s\n" % "\t".join(newd));
        