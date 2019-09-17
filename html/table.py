import sys
import os
import dominate
from dominate.tags import *
from pybedtools import BedTool
from dominate.util import raw

#session = "759704951_gifQr7fJt2QL5GNdnr0jGn5Mz7AP"
#session = "759704951"

def add_ucsc(interval, session):
    ucsc_link = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=%s%%3A%d%%2D%d&hgsid=%s' % (interval.chrom, interval.start, interval.stop, session)
    return ucsc_link
    #return "\n<a href=%s>\n</a>\n" % ucsc_link




peaks = BedTool(sys.argv[1])
headers = ['ucsc', 'chrom', 'start', 'end', 'strand', 'gene', 'score', 'max coverage']
#dtypes = ['string', 'numeric', 'numeric', 'string', 'string', 'numeric', 'numeric']
dtypes = [0, 0, 1, 1, 0, 0, 1, 1]

#print(add_ucsc(peaks[0]));
#sys.exit()



doc = dominate.document(title='test table')

with open(sys.argv[2]) as f:
    _style = f.read()
    
with open(sys.argv[3]) as f:
    plain_script = f.read()



with doc.head:
    style(_style)
    #link(rel='stylesheet', href=sys.argv[2])
    #script(type='text/javascript', src=sys.argv[3])
    #_style = style()
    #_style.add(table(border_spacing=0, width=100, border="1px solid"))

with doc:
    p(strong("lets try "))
    p("to learn a bit")
    input(type="text", id="myInput", onkeyup="my_search(4)", placeholder="Search for genes..")
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x[1][0], onclick='sortTable(%d, %d)' % (x[0], x[1][1])) for x in enumerate(zip(headers, dtypes))  ])
        for interval in peaks:
            with tr():
                td(raw('<a href=%s target="_blank">ucsc_link</a>' % add_ucsc(interval, sys.argv[4])) )
                td(interval.chrom)
                td(interval.start)
                td(interval.stop)
                td(interval.strand)
                td(interval.attrs['start_gene'])
                td(interval.score)
                td(interval.attrs['topcoverage'])
                
    _script = script(type='text/javascript')
    _script.add_raw_string(plain_script)
    
        

    
print(doc.render());

    
    

