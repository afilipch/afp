import sys
import os
import dominate
import argparse
from dominate.tags import *
from pybedtools import BedTool
from dominate.util import raw

from afbio.html.methods import add_ucsc

parser = argparse.ArgumentParser(description='Creates html table for the annotated ChIP binding peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated ChIP peaks, gff format");
parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
parser.add_argument('--name', nargs = '?', default="unknown", type = str, help = "Name of the sample");
parser.add_argument('--top', nargs = '?', default=300, type = int, help = "Shows only [top] peaks");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory (tsv and html files will be written there)")
args = parser.parse_args();

chr_dict = {'NC_003450.3': 'chr1', 'pJC1-Plys::GntR': 'chr2'}


def flatten_peak(peak):
    return peak.chrom, peak.start, peak.stop, peak.name, peak.strand, peak.attrs["gene"], peak.attrs["genesymbol"], peak.attrs["gtype"], peak.attrs["tss"], peak.attrs["atg"], peak.attrs["annotation"], peak.attrs["function"], peak.attrs["topcoverage"], peak.attrs["area_coverage"]
    

_title = "Binding peaks for the sample: %s" % args.name
peaks = list(sorted(BedTool(args.path), key = lambda x: float(x.score), reverse = True))[:args.top]

header_main = ["ucsc", 'chrom', 'start', 'end', 'top', 'strand']
dtypes_main = [0, 0, 1, 1, 1, 0]
header_attrs =  ['geneID', 'genesymbol', 'gtype', 'tss', 'atg', 'annotation', 'function', 'topcoverage', 'area_coverage']
dtypes_attrs = [0, 0, 0, 1, 1, 0, 0, 1, 1]
header = header_main + header_attrs
dtypes = dtypes_main + dtypes_attrs





doc = dominate.document(title=_title)

with open(args.css) as f:
    _style = f.read()
    
with open(args.js) as f:
    plain_script = f.read()



with doc.head:
    style(_style)


with doc:
    p(strong(_title))
    input(type="text", id="inp1", onkeyup='my_search(7, "inp1", "myTable")', placeholder="Search for a genesymbol..")
    input(type="text", id="inp2", onkeyup='my_search(6, "inp2", "myTable")', placeholder="Search for a geneID..")
    input(type="text", id="inp3", onkeyup='my_search(8, "inp3", "myTable")', placeholder="Filter for target type..")
    input(type="number", id="inp4", onkeyup='my_filter_greater(13, "inp4", "myTable")', placeholder="Filter top coverage content greater than..")
    input(type="number", id="inp5", onkeyup='my_filter_lesser_abs(9, "inp5", "myTable")', placeholder="Filter distance to TSS closer than..")
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x[1][0], onclick='sortTable(%d, %d)' % (x[0], x[1][1])) for x in enumerate(zip(header, dtypes))  ])
        for peak in peaks:
            with tr():
                td(raw('<a href=%s target="_blank">ucsc_link</a>' % add_ucsc(peak, args.ucsc, flank=25, chr_dict=chr_dict)) )
                for entry in flatten_peak(peak):
                    td(entry)
                #td(peak.chrom)
                #td(peak.start)
                #td(peak.stop)
                #td(peak.name)
                #td(peak.strand)
                #td(peak.attrs["gene"])
                #td(peak.attrs["genesymbol"])
                #td(peak.attrs["gtype"])
                #td(peak.attrs["tss"])
                #td(peak.attrs["atg"])
                #td(peak.attrs["annotation"])
                #td(peak.attrs["function"])
                #td(peak.attrs["topcoverage"])
                #td(peak.attrs["area_coverage"])
                
                
    _script = script(raw(plain_script), type='text/javascript')
    
        

with open(os.path.join(args.outdir, 'peaks.html'), 'w') as f:
    #print(doc.render())
    f.write(doc.render());
    

with open(os.path.join(args.outdir, 'peaks.tsv'), 'w') as f:
    f.write("%s\n" % "\t".join(header[1:]));
    for peak in peaks:
        f.write("%s\n" % "\t".join([str(x) for x in flatten_peak(peak)]));
        

    


    
    

