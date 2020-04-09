import sys
import os
from collections import defaultdict
import dominate
import argparse
from dominate.tags import *
from dominate.util import raw
from pybedtools import BedTool


from afbio.html.methods import add_ucsc

parser = argparse.ArgumentParser(description='Creates html table for the clustered genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the clustered differential table, tsv format");
parser.add_argument('--annotation', nargs = '?', required=True, type = str, help = "Path to the gene annotation");
parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
parser.add_argument('--num_conditions', nargs = '?', required = True, type = int, help = "Number of conditions");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory (tsv and html files will be written there)")
parser.add_argument('--genes_dir', nargs = '?', required=True, type = str, help = "Path to the directory with genes html reports")
parser.add_argument('--plot_dir', nargs = '?', required=True, type = str, help = "Path to the directory with trend plots")
args = parser.parse_args();






    
def annotate(d, gene2annotation):
    gene = gene2annotation[d[0]]
    name = gene.attrs['genesymbol']
    if(name == 'None'):
        name = gene.attrs['ID']
    ann = gene.attrs['product']
    return [name, ann] + d[1:], gene


gene2annotation = dict([ (x.attrs['ID'], x) for x in BedTool(args.annotation)])

data_dict = defaultdict(list)
with open(args.path) as f:
    header = next(f).strip().split("\t")
    cl_index = header.index('cluster')
    header = header[:cl_index+args.num_conditions+1]
    header.insert(1, 'UCSC')
    header.insert(2, 'product')
    
    for l in f:
        a = l.strip().split("\t")[:cl_index+args.num_conditions+1]
        data_dict[int(a[cl_index])].append(a)
        

with open(args.css) as f:
    _style = f.read()
    
with open(args.js) as f:
    plain_script = f.read()



for cl in sorted(data_dict.keys()):
    data = data_dict[cl]
    data.sort(key = lambda x: float(x[2]), reverse=True)
    
    _title = "Gene expression trends: cluster number %d" % cl
    doc = dominate.document(title=_title)
    with doc.head:
        style(_style)
    with doc:
        p(strong(_title))
        with div(cls="row", style="display: flex"):
            img(src = os.path.join(args.plot_dir, 'trend_c%d.png' % cl), alt="trend", style="width:75%")
        br()
        br()
        input(type="text", id="inp1", onkeyup='my_search(1, "inp1", "myTable")', placeholder="Search for a genesymbol..")
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
    
    path = os.path.join(args.outdir, "cluster_%d" % cl)
    with open(path + ".html" , 'w') as f:
        f.write(doc.render());
        
    #with open(path + ".tsv" , 'w') as f:
        #f.write("%s\n" % "\t".join(header[1:]));
        #for d in data:
            #newd, gene = annotate(d, gene2annotation)
            #f.write("%s\n" % "\t".join(newd));
        
