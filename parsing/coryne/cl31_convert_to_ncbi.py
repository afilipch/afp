#! /usr/bin/python
'''Converts transcripts raw file into ncbi format'''

import argparse
import sys
import os
from collections import defaultdict, Counter

from pybedtools import BedTool
from afbio.pybedtools_af import construct_gff_interval



parser = argparse.ArgumentParser(description='Converts transcripts raw file into ncbi format');
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts, raw format");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format");
args = parser.parse_args();

with open(args.genome) as f:
    chrom = next(f).strip()[1:]

with open(args.transcripts) as f:
    next(f);
    for l in f:
        a = [x.strip().replace('"', '').replace(';', '') for x in l.strip().split("\t")]
        aint = construct_gff_interval( chrom, int(a[2]), int(a[3]), 'gene', score='0', strand=a[4], source='ncbi_af', frame='.', attrs=[ ('ID', 'gene-%s' % a[0]), ('gene_biotype', 'protein_coding'), ('Name', 'gene-%s' % a[0]), ('product', a[1]) ]  );
        sys.stdout.write(str(aint))
        
#453-Cg-phage-CL31_S7_L001_R1_001_contig_1       un      gene    2       2422    0       +       .       ID=gene-CL_1_Hyp.; gene_biotype=protein_coding; Name=gene-CL_1_Hyp.; product=Protein (A. faecalis)
# ['CL_53', 'Hyp. Protein (Corynebacterium)', '44425', '44796', '-']


#aint = construct_gff_interval(parent.chrom, parent.start, parent.stop, 'gene', score='0', strand=parent.strand, source='ncbi_af', frame='.', attrs=[('ID', parent.attrs['ID']), ('genesymbol', parent.attrs.get('gene', 'None')), ('annotation', ga.attrs.get('Note', 'None')), ('product', ga.attrs.get('product', 'None')), ('cds', "%d:%d" % (ga.start+1, ga.stop)), ('tss_variants', '1'), ('distance', str(distance))  ] );
    
    
    
    
    
    
    
    



