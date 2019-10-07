import sys
import os
import dominate
import argparse
from dominate.tags import *
from pybedtools import BedTool
from dominate.util import raw
from Bio import SeqIO

from afbio.html import add_ucsc
from math import log

parser = argparse.ArgumentParser(description='Filters the detected peaks based on the distribution of their scores');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the AT-rich areas, gff format");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format");
parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
parser.add_argument('--flank', nargs = '?', default=60, type = int, help = "Peak plank length");
args = parser.parse_args();

PAM_LENGTH = 21;
HALF_PAM = 10;


def split_sequence(peak, genome, flank):
    seq = peak2sequence(peak, genome, flank)

    sense, antisense = [seq], [seq]
    
    sense_start = flank + int(peak.attrs['pam_sense']) - HALF_PAM - 3
    sense_stop = sense_start + PAM_LENGTH + 2
    if(sense_start>=0 and sense_stop<=flank*2):
        sense = [seq[:sense_start], seq[sense_start:sense_start+3], seq[sense_start+3:sense_stop], seq[sense_stop:]]
    
    antisense_start = flank - int(peak.attrs['pam_antisense']) - HALF_PAM 
    antisense_stop = antisense_start + PAM_LENGTH + 2
    if(antisense_start>=0 and antisense_stop<=flank*2):
        antisense = [seq[:antisense_start], seq[antisense_start:antisense_stop-3], seq[antisense_stop-3:antisense_stop], seq[antisense_stop:] ] 
     
    cl = flank//4
    tip = ["-"*(flank-cl), '*'*(cl*2), "-"*(flank-cl)]
    
    sense[-1] += '\n'
    tip[-1] += '\n'
    #sys.stderr.write(str(tip))
    
    return sense, antisense, tip
    

def peak2sequence(peak, genome, flank):
    center = int(peak.name)
    start = max(center - flank, 0)
    stop = center + flank
    
    if(peak.strand == '+'):
        return str(genome[peak.chrom][start:stop].seq.upper())
    elif(peak.strand == '-'):
        return str(genome[peak.chrom][start:stop].seq.reverse_complement().upper()) 


def set_pam_score(peak):
    tss = 1+abs(int(peak.attrs['tss']))
    if(tss<60):
        tss = 10;
    zscore = float(peak.score)
    pam = int(peak.attrs['pam_min'])
    if(pam<=11):
        pam = 5;
    return (1/tss)*zscore*(1/pam**2)*1000



peaks = BedTool(args.path)
genome = SeqIO.to_dict(SeqIO.parse(args.genome, 'fasta'))
newpeaks = []
for peak in peaks:
    score = set_pam_score(peak)
    peak.attrs['pam_score'] = "%1.1f" % score
    newpeaks.append(peak);
    #split_sequence(peak, genome, args.flank)


peaks = sorted(newpeaks, key=lambda x: float(x.attrs['pam_score']), reverse=True)




cgps_dict = {'in': 'in', 'unk': 'close', 'out': 'out'}
headers = ['ucsc', 'start', 'end', 'strand', 'gene', 'TSS distance', 'peak score', 'pam distance', 'pam score', 'sequence 5->3']
dtypes = [0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0]
sense_colors = ['black', 'red', 'purple', 'black']
antisense_colors = ['black', 'purple', 'red', 'black']
tip_colors = ['gray', 'black', 'gray']


doc = dominate.document(title='CgpS peaks candidates for dCas9 countersilencing')

with open(args.css) as f:
    _style = f.read()
    
with open(args.js) as f:
    plain_script = f.read()



with doc.head:
    style(_style)


with doc:
    p(strong("CgpS peaks candidates for dCas9 countersilencing"))
    input(type="text", id="myInput", onkeyup="my_search(5)", placeholder="Search for a gene ...")
    #input(type="number", id="myInputGreater", onkeyup="my_filter_greater(6)", placeholder="Filter AT content greater than..")
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x[1][0], onclick='sortTable(%d, %d)' % (x[0], x[1][1])) for x in enumerate(zip(headers, dtypes))  ])
        for interval in peaks:
            center = int(interval.name)
            sense, antisense, tip = split_sequence(interval, genome, args.flank)
            #print(sense)
            with tr():
                td(raw('<a href=%s target="_blank">ucsc_link</a>' % add_ucsc(interval, args.ucsc)) )
                td(center-args.flank)
                td(center+args.flank)
                td(interval.strand)
                td(interval.attrs['gene'])
                td(interval.attrs['tss'])
                td(interval.score)
                td(interval.attrs['pam_min'])
                td(interval.attrs['pam_score'])
                with td(__pretty=False, style= 'font-family:monospace; font-size: 12pt'):
                    for seq, color in zip(sense, sense_colors):
                        span(seq, style="color:%s" % color)
                    for seq, color in zip(tip, tip_colors):
                        span(seq, style="color:%s" % color)
                    for seq, color in zip(antisense, antisense_colors):
                        span(seq, style="color:%s" % color)
                #td(antisense[0]);
                
    _script = script(type='text/javascript')
    _script.add_raw_string(plain_script)
    
        

    
print(doc.render());

    
    

