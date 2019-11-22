#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Annotates the discovered peaks'''

import argparse
import sys
import os
from collections import defaultdict
from bisect import bisect_right, bisect_left
from afbio.sequencetools import coverage2dict


import pandas as pd;
import numpy as np;
import matplotlib.pyplot as plt;
from pybedtools import BedTool

from afbio.pybedtools_af import construct_gff_interval

parser = argparse.ArgumentParser(description='Annotates the discovered peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the detected peaks");
parser.add_argument('--genes', nargs = '?', required=True, type = str, help = "Path to the gene annotation file");
parser.add_argument('--coverage', nargs = '?', type = str, help = "Path to the coverage track, bed format");
parser.add_argument('--overlap', nargs = '?', default=0.5, type = float, help = "Min overlap required overlap between a peak and genomic feature (as fraction of peak length)");
parser.add_argument('--maxshift', nargs = '?', default=1, type = int, help = "Max allowed shift (in nucleotides) of the peak top position downstream to start of the gene, to be still counted as peak upstream the gene");
parser.add_argument('--promoter_distance', nargs = '?', default=700, type = int, help = "Max allowed distance to the gene start for the peak to be reported as promoter");
parser.add_argument('--flen', nargs = '?', default=50, type = int, help = "Length of the peak\'s flanks to be included into analyses");
parser.add_argument('--custom', nargs = '?', default=False, const=True, type = bool, help = "If set the annotation genes are supposed to be already processed, if not they are supposed to be in NCBI gff3 format");
args = parser.parse_args();


if(os.stat(args.path).st_size != 0):


    #################################################################################################################################################################
    ### Read the input

    genes = BedTool(args.genes);
    peaks = BedTool(args.path);
    offset = len(peaks[0].fields)

    #################################################################################################################################################################
    ### Convert NCBI genes the input
    if(not args.custom):
        coding_genes = [];
        
        def genebox2coding(genebox):
            name = genebox[0].name
            function = genebox[1].attrs.get('Note', 'None').replace(';', ':');
            annotation = genebox[1].attrs['product'].replace(';', ':')
            return construct_gff_interval('chr1', genebox[0].start, genebox[0].end, 'protein_coding', score='0', strand=genebox[0].strand, source='.', frame='.', attrs=[ ('Name', name), ('function', function), ('annotation', annotation)])
        
        genebox = [];
        for interval in genes:
            if(interval[2] in ['gene', 'pseudogene']):
                if(genebox):
                    if(len(genebox)==2):
                        coding_genes.append(genebox2coding(genebox))
                        #print([x[2] for x in genebox])
                genebox = [interval]
            elif(genebox):
                genebox.append(interval)
        else:
            if(len(genebox)==2):
                coding_genes.append(genebox2coding(genebox))
                
    #for interval in coding_genes:
        #print(interval)
        genes = coding_genes;
        

    genes2annotation = dict([ (x.name, (x.attrs['annotation'], x.attrs['function']) ) for x in genes])
                            

    #################################################################################################################################################################
    ### Get names of the overlapping genes

    
    def get_gene_name(intersection, offset):
        #print([x.strip().split('=') for x in intersection[offset+8].strip(';').split(";")])
        attrs = dict( [x.strip().split('=') for x in intersection[offset+8].strip(';').split(";")])
        return attrs['Name']

    peak2genenames = defaultdict(list);
    for el in peaks.intersect(genes, wo = True, f = args.overlap):
        peak2genenames[el.name].append(get_gene_name(el, offset))
        
        
    if(peaks.file_type == 'gff'):
        temp_peaks = []
        for interval in peaks:
            interval.attrs['genes'] = ",".join(peak2genenames.get(interval.name, ['None']))
            temp_peaks.append(interval)
        peaks = temp_peaks;
    else:
        temp_peaks = []
        for interval in peaks:
            top = int(interval.name)
            start = max((top-args.flen,0))
            end = top+args.flen +1
            anint = construct_gff_interval(interval.chrom, start, end, 'peak', score=interval.score, strand=interval.strand, source='af_peak_detection', frame='.', attrs=[ ('Name', interval.name), ('genes', ",".join(peak2genenames.get(interval.name, ['None'])))])
            temp_peaks.append(anint)
        peaks = temp_peaks;
            
        
        

        
    #################################################################################################################################################################
    ### Get coverage annotation
    
    if(args.coverage):
        coverage = coverage2dict(args.coverage)
        cov_list = [coverage2dict(args.coverage, cpos=x) for x in range(3, 7)]
        for interval in peaks:
            interval.attrs['topcoverage'] =  "%1.3f" % coverage[interval.chrom][int(interval.name)]
            interval.attrs['other_coverage'] =  ",".join(["%1.3f" % x[interval.chrom][int(interval.name)] for x in cov_list])
            


    
    #sys.exit()
        

    
    #################################################################################################################################################################    
    ###Get closest genomic starts upstream to the dected peaks
    #plusstarts = [x.start for x in genes if gene.strand == '+']
    #minusstarts = [x.end-1 for x in genes if gene.strand == '-']

    genestarts = [ (x.start, x.strand, x.name) if x.strand == '+' else (x.end-1, x.strand, x.name) for x in genes]


    def findclosest(interval, genestarts, maxshift):
        pos = int(interval.name)
        raw_distances = [x[0]-pos if x[1] == '+' else pos-x[0] for x in genestarts];
        distances = [abs(x) if x > -maxshift else 10**8 for x in raw_distances]
        minindex = np.argmin(distances);
        start, strand, name = genestarts[minindex];
        return distances[minindex], name, strand
    
    
    def add_type(interval, maxdistance):
        if(interval.attrs['genes'].split(",")[0] != 'None'):
            return "Gene"
        elif(int(interval.attrs['start_gene_distance']) <= maxdistance):
            return "Promoter"
        else:
            return "Intergenic"
        

                
        
    peak2genestarts = {}
    for interval in peaks:
        ss_distance, ss_genename, ss_strand = findclosest(interval, genestarts, args.maxshift)
        interval.attrs['start_gene'] = ss_genename
        interval.attrs['start_gene_distance'] =  "%d" % ss_distance
        interval.attrs['start_gene_strand'] = ss_strand
        interval.attrs['start_annotation'], interval.attrs['start_function'] = genes2annotation[ss_genename]
        interval.attrs['type'] = add_type(interval, args.promoter_distance)
        sys.stdout.write(str(interval))







