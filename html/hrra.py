import sys
import os
import dominate
import argparse
from dominate.tags import *
from pybedtools import BedTool
from dominate.util import raw

parser = argparse.ArgumentParser(description='Creates html table for HrrA all-in table');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "HrrA all-in table, tsv format");
parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
args = parser.parse_args();

#https://genome.ucsc.edu/s/philipchick/cgps?position=chr10:69,644,222-69,644,999



def add_ucsc(start, stop, session):
    ucsc_link = 'https://genome.ucsc.edu/s/%s?position=chr1:%s-%s' % (session, start, stop)
    return ucsc_link



all_in_table = []
with open(args.path) as f:
    headers = next(f).strip().split("\t")
    for l in f:
        l = l.strip()
        if(l):
            all_in_table.append(l.split("\t"))
        
#print(headers)
#sys.exit()

headers = ["ucsc", 'start', 'end', 'Gene ID', 'Gene symbol', 'Distance ATG', 'Distance to TSS', 'ChAP T=pre', 'ChAP T=0h', 'ChAP T=30m', 'ChAP T=2h', 'ChAP T=4h', 'ChAP T=9h', 'ChAP T=24h', 'mRNA wt T=0h', 'mRNA wt T=30m', 'mRNA wt T=4h', 'mRNA DhrrA T=0h', 'mRNA DhrrA T=30m', 'mRNA DhrrA T=4h', 'Log2 DhrrA/WT T=0h', 'Log2 DhrrA/WT T=30m', 'Log2 DhrrA/WT T=4h', 'Predicted Function', 'Annotation', "Divergent"]

dtypes = [0, 1, 1, 0, 0] + [1]*(len(headers)-7) + [0, 0, 0]
#sys.stderr.write(str(dtypes))
ndict = {'None': 0}
#print(add_ucsc(peaks[0]));
#sys.exit()



doc = dominate.document(title='HrrA All-In table')

with open(args.css) as f:
    _style = f.read()
    
with open(args.js) as f:
    plain_script = f.read()



with doc.head:
    style(_style)


with doc:
    p(strong("HrrA All-In table"))
    input(type="text", id="myInput", onkeyup="my_search(4)", placeholder="Search for a genesymbol..")
    #input(type="number", id="myInputGreater", onkeyup="my_filter_greater(6)", placeholder="Filter AT content greater than..")
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x[1][0], onclick='sortTable(%d, %d)' % (x[0], x[1][1])) for x in enumerate(zip(headers, dtypes))  ])
        for a in all_in_table:
            with tr():
                #td(raw('<a href=%s target="_blank">ucsc_link</a>' % add_ucsc(interval, args.ucsc)) )
                td(raw('<a href=%s target="_blank">ucsc_link</a>' % add_ucsc(a[1], a[2], args.ucsc)) )
                for el in a[1:]:
                    td(ndict.get(el, el))
                
    _script = script(type='text/javascript')
    _script.add_raw_string(plain_script)
    
        

    
print(doc.render());

    
    

