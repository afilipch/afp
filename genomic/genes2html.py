import sys
import os
import argparse
from collections import defaultdict
from itertools import groupby

from pybedtools import BedTool;

import dominate
from dominate.tags import *
from dominate.util import raw
from Bio import SeqIO

from afbio.html.methods import add_ucsc
from math import log
from afbio.pybedtools_af import interval2seq

parser = argparse.ArgumentParser(description='Generates html pages for the provided annotated genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated genes");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");
parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
args = parser.parse_args();

### SET CONSTANTS ###

PHAGE = {'0': 'No', '1': 'Yes', None: 'Unknown'}


### FUNCTIONS ###

def keyfunc(transcript):
    return transcript.name

def transcripts2gene(group):
    first = group[0]
    a = [int(x) for x in first.attrs['cds'].split(":")]
    first.start, first.stop = a[0]-1, a[1]
    
    starts = sorted(list(set([x.start for x in group])))
    stops = sorted(list(set([x.stop for x in group])))
    if(first.strand == '-'):
        starts, stops = stops, starts
        
    if(len(starts)>1):
        print(first.name)
        
        
    first.attrs["alt_tss"] = ",".join([str(x) for x in starts]);
    first.attrs["alt_tts"] = ",".join([str(x) for x in stops]);
    
    return first 

    



transcripts = list(BedTool(args.path))
transcripts.sort(key = keyfunc);
groups = groupby(transcripts, keyfunc);

genes = []
for key, data in groups:
    genes.append(transcripts2gene(list(data)));
    
genome = SeqIO.to_dict(SeqIO.parse(args.genome, 'fasta'))
#interval2seq(interval, reference)
#print(genes[0]);
#sys.exit()
#for gene in genes:
    #print(gene)
    #print()
    #for kv in gene.attrs.items():
        #print(kv)
    #break;



with open(args.css) as css:
    _style = css.read()

for gene in genes:
########### HTML SECTION ###########
    with open(os.path.join(args.outdir, "%s.html" % gene.name), 'w') as f:
        doc = dominate.document(title="%s: gene overview" % gene.name)
        with doc.head:
            style(_style)
        with doc:
            p(strong(gene.name))
            br()
            with p():
                b("Gene Symbol:")
                raw("\t%s" % gene.attrs['genesymbol'])
            with p():
                b("CDS:")
                raw('<a href=%s target="_blank">\t%s(%s) %d-%d</a>' % (add_ucsc(gene, args.ucsc), gene.chrom, gene.strand, gene.start, gene.stop) )
                #raw("\t%s(%s) %d-%d " % (gene.chrom, gene.strand, gene.start, gene.stop))
            with p():
                b("Alternative TSS:")
                raw("\t%s" % gene.attrs["alt_tss"].replace(",", ", ") )
            with p():
                b("Annotation:")
                raw("\t%s" % gene.attrs['annotation'] )
            with p():
                b("Function:")
                raw("\t%s" % gene.attrs.get('function', gene.attrs.get('product', 'None')) )
            with p():
                b("Phage:")
                raw("\t%s" % PHAGE[gene.attrs.get('phage')] )
            with p():
                b("Sequence:")
                raw("\t%s" %  interval2seq(gene, genome))
                

                    
                
                


        f.write(doc.render());
    #break;
