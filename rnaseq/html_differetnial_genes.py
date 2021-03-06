import sys
import os
import dominate
import argparse
from dominate.tags import *
from dominate.util import raw
from pybedtools import BedTool

from afbio.html.methods import add_ucsc, ucsc_convert_chromosomes

parser = argparse.ArgumentParser(description='Creates html table for the differentially expressed genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the assigned differential table, tsv format");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the gene annotation");
parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
parser.add_argument('--name', nargs = '?', default="unknown", type = str, help = "Name of the experiment");
parser.add_argument('--top', nargs = '?', default=300, type = int, help = "Shows only [top] changed genes");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory (tsv and html files will be written there)")
parser.add_argument('--genes_dir', nargs = '?', required=True, type = str, help = "Path to the directory with genes html reports")
parser.add_argument('--chrom_dict', nargs = '+', type = str, help = "UCSC chrom dictionary. Format real_chrom,ucsc_chrom")
args = parser.parse_args();





    

_title = "Differential gene expression of: %s" % args.name

data = []
with open(args.path) as f:
    header = next(f).strip().split("\t")
    for l in f:
        data.append(l.strip().split("\t"))
        
data.sort(key = lambda x: float(x[2]), reverse=True)

data = data[:args.top]

annotation = BedTool(args.annotation)
if(args.chrom_dict):
    chrom_dict = dict([ tuple(x.split(",")) for x in args.chrom_dict ])
else:
    chrom_dict = ucsc_convert_chromosomes(annotation)
#sys.stderr.write('%s\n' % chrom_dict)
gene2annotation = dict([ (x.attrs['ID'], x) for x in annotation])


def annotate(d, gene2annotation):
    name = d[0];
    gene = gene2annotation[d[0]]
    
    alt_name = gene.attrs.get('cg', 'none')
    genesymbol = gene.attrs.get('genesymbol', name)
    product = gene.attrs['product']
    return [name, alt_name, genesymbol, product] + d[1:], gene
    






doc = dominate.document(title=_title)

with open(args.css) as f:
    _style = f.read()
    
with open(args.js) as f:
    plain_script = f.read()



with doc.head:
    style(_style)

header = header[:1] + ['alt_name', 'genesymbol', 'UCSC', 'product'] + header[1:]



with doc:
    p(strong(_title))
    input(type="text", id="inp1", onkeyup='my_search(0, "inp1", "myTable")', placeholder="Search for a genesymbol..")
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x) for x in header])
        for d in data:
            newd, gene = annotate(d, gene2annotation)
            with tr():
                td(raw('<a href=%s target="_blank">%s</a>' % (os.path.join(args.genes_dir, newd[0] + ".html"), newd[0]) ))
                td(newd[1]);
                td(newd[2]);
                td(raw('<a href=%s target="_blank">link</a>' % add_ucsc(gene, args.ucsc, chr_dict=chrom_dict) ))
                for el in newd[3:]:
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
    f.write("%s\n" % "\t".join(header[:3] + header[4:]));
    for d in data:
        newd, gene = annotate(d, gene2annotation)
        f.write("%s\n" % "\t".join(newd));
        
