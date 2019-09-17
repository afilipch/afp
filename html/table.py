import sys
import os
import dominate
from dominate.tags import *
from pybedtools import BedTool


peaks = BedTool(sys.argv[1])
headers = ['chrom', 'start', 'end', 'strand', 'gene', 'score', 'max coverage']
#dtypes = ['string', 'numeric', 'numeric', 'string', 'string', 'numeric', 'numeric']
dtypes = [0, 1, 1, 0, 0, 1, 1]


doc = dominate.document(title='test table')

with doc.head:
    link(rel='stylesheet', href=sys.argv[2])
    script(type='text/javascript', src=sys.argv[3])
    #_style = style()
    #_style.add(table(border_spacing=0, width=100, border="1px solid"))

with doc:
    p(strong("lets try "))
    p("to learn a bit")
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x[1][0], onclick='sortTable(%d, %d)' % (x[0], x[1][1])) for x in enumerate(zip(headers, dtypes))  ])
        for interval in peaks:
            with tr():
                td(interval.chrom)
                td(interval.start)
                td(interval.stop)
                td(interval.strand)
                td(interval.attrs['start_gene'])
                td(interval.score)
                td(interval.attrs['topcoverage'])
                
                
        

    
print(doc);

    
    

