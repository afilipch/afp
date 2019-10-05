import sys
import os
import dominate
import argparse
from dominate.tags import *
from pybedtools import BedTool
from dominate.util import raw

parser = argparse.ArgumentParser(description='Filters the detected peaks based on the distribution of their scores');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the AT-rich areas, gff format");
parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
args = parser.parse_args();

#https://genome.ucsc.edu/s/philipchick/cgps?position=chr10:69,644,222-69,644,999



def add_ucsc(interval, session):
    ucsc_link = 'https://genome.ucsc.edu/s/%s?position=chr1:%d-%d' % (session, interval.start-25, interval.stop+25)
    return ucsc_link





peaks = BedTool(args.path)
cgps_dict = {'in': 'in', 'unk': 'close', 'out': 'out'}
headers = ['ucsc', 'chrom', 'start', 'end', 'CgpS', 'upstream', 'max AT content [%]', 'flank AT content [%]']
#dtypes = ['string', 'numeric', 'numeric', 'string', 'string', 'numeric', 'numeric']
dtypes = [0, 0, 1, 1, 0, 0, 1, 1]

#print(add_ucsc(peaks[0]));
#sys.exit()



doc = dominate.document(title='AT rich areas')

with open(args.css) as f:
    _style = f.read()
    
with open(args.js) as f:
    plain_script = f.read()



with doc.head:
    style(_style)


with doc:
    p(strong("AT rich areas"))
    input(type="text", id="myInput", onkeyup="my_search(4)", placeholder="Filter relation to CgpS..")
    input(type="number", id="myInputGreater", onkeyup="my_filter_greater(6)", placeholder="Filter AT content greater than..")
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x[1][0], onclick='sortTable(%d, %d)' % (x[0], x[1][1])) for x in enumerate(zip(headers, dtypes))  ])
        for interval in peaks:
            with tr():
                td(raw('<a href=%s target="_blank">ucsc_link</a>' % add_ucsc(interval, args.ucsc)) )
                td(interval.chrom)
                td(interval.start-25)
                td(interval.stop+25)
                td(cgps_dict[interval.attrs['type']])
                td(interval.attrs['upstream'])
                td("%1.1f" % (float(interval.score)*100))
                td("%1.1f" % (float(interval.attrs['at_flank'])*100))
                
    _script = script(type='text/javascript')
    _script.add_raw_string(plain_script)
    
        

    
print(doc.render());

    
    

