import sys
import os
import dominate
import argparse
from collections import defaultdict
from dominate.tags import *
from dominate.util import raw


from afbio.generators import get_only_files

parser = argparse.ArgumentParser(description='Generates html report of CHAP analyses by chipchap.py');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the differential folder");
parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--name', nargs = '?', default="experiment", type = str, help = "Name of the experiment");
parser.add_argument('--tables', nargs = '?', required=True, type = str, help = "path to the folder with separate (per pair) html tables");
args = parser.parse_args();

### SET CONSTANTS ###

labels = ["Condition 1", "Condition 2", "Num of greater peaks", "Link"] #, "Mapped non-uniquely", "Unmapped reads", "Mapped Discordantly"]
dtypes = [0,0,1,0]

files = get_only_files(args.path)


data = [];
for f in files:
    c1, c2 = os.path.basename(f).split(".")[0].split("_GREATER_")
    count = len(open(f).readlines())
    ref = os.path.join(args.tables, os.path.basename(f).split(".")[0]) + ".html"
    data.append((c1,c2,count, ref))
    

#references = [os.path.join(args.tables, os.path.basename(x).split(".")[0]) + ".html" for x in files]




        





########### HTML SECTION ###########
doc = dominate.document(title="%s: differential CHAP analyses overview" % args.name)
with open(args.css) as f:
    _style = f.read()
with doc.head:
    style(_style)
    
with open(args.js) as f:
    plain_script = f.read()
    
with doc:
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x[1][0], onclick='sortTable(%d, %d)' % (x[0], x[1][1])) for x in enumerate(zip(labels, dtypes))  ])
        for d in data:
            with tr():
                td(d[0]) 
                td(d[1])
                td(d[2])
                td(raw('<a href=%s target="_blank">%s</a>' % (d[3], "click me") )) 



                
    _script = script(raw(plain_script), type='text/javascript')

print(doc.render());





