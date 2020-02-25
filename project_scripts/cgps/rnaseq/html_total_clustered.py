import sys
import os
from collections import defaultdict
import dominate
import argparse
from dominate.tags import *
from dominate.util import raw
from pybedtools import BedTool


from afbio.html.methods import add_ucsc

parser = argparse.ArgumentParser(description='Creates global html table for the clusters');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the clustered differential table, tsv format");
parser.add_argument('--html_dir', nargs = '?', required=True, type = str, help = "Path to the directory with clusters html reports")
parser.add_argument('--plot_dir', nargs = '?', required=True, type = str, help = "Path to the directory with trend plots")
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
args = parser.parse_args();






    


data_dict = defaultdict(int)
with open(args.path) as f:
    header = next(f).strip().split("\t")
    cl_index = header.index('cluster')
    for l in f:
        a = l.strip().split("\t")
        data_dict[int(a[cl_index])] += 1
        

with open(args.css) as f:
    _style = f.read()
    





header = "Cluster ID", "Number of genes", "Trend"
_title = "Gene expression trends"
doc = dominate.document(title=_title)
with doc.head:
    style(_style)
with doc:
    p(strong(_title))
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x) for x in header])
        for cl in sorted(data_dict.keys()):
            count = data_dict[cl]
            with tr():
                td(raw('<a href=%s target="_blank">%s</a>' % (os.path.join(args.html_dir, "cluster_%d.html" % cl), cl) ))
                td(count)
                td(img(src = os.path.join(args.plot_dir, 'trend_c%d.png' % cl), alt="trend", style="width:40%"))
print(doc.render());
        
    #with open(path + ".tsv" , 'w') as f:
        #f.write("%s\n" % "\t".join(header[1:]));
        #for d in data:
            #newd, gene = annotate(d, gene2annotation)
            #f.write("%s\n" % "\t".join(newd));
        
