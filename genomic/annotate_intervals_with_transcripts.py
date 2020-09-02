#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Annotates genomic intervals with transcripts'''

import argparse
import sys
import os
from collections import defaultdict, Counter
from afbio.sequencetools import coverage2dict


import numpy as np;
import matplotlib.pyplot as plt;
from pybedtools import BedTool

from afbio.pybedtools_af import construct_gff_interval, read_comments


parser = argparse.ArgumentParser(description='Annotates genomic intervals with transcripts');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genomic regions, gff format");
parser.add_argument('--transcripts', required = True, nargs = '?',  type = str, help = "Path to the transcripts regions, gff format");
parser.add_argument('--inside', nargs = '?', default=200, type = int, help = "Maximum allowed distance to TSS while inside a gene");
parser.add_argument('--maxd', nargs = '?', default=800, type = int, help = "Maximum allowed distance to TSS");
parser.add_argument('--format', nargs = '?', default='svg', type = str, help = "Plot format, svg by default");
parser.add_argument('--outdir', nargs = '?', required=False, type = str, help = "Path to the output plot directory")
args = parser.parse_args();


if(os.stat(args.path).st_size == 0):
    sys.exit("###annotate\nInput file is empty, empty output is produced\n")

STUB_TR = construct_gff_interval( "unknown", 0, 2, 'fake', score='0', strand="+", source='ff', frame='.', attrs= [("Name", "fake"), ("annotation", "None"), ("function", "None"), ("genesymbol", "fake"), ("distance", "NaN")])

def get_anti(tr, mind):
    return [('anti_gene', tr.name), ('anti_genesymbol', tr.attrs['genesymbol']),  ('anti_tss', str(mind))]
    

def annotate_position(interval, tr_dict, maxd, inside):
    center = int((interval.stop + interval.stop)/2)
    tr_plus, tr_minus = tr_dict.get(interval.chrom, [None, None])
    if(not tr_plus):
        mindistance = float("NaN")
        gtype = "intergenic"
        transcript = STUB_TR
        anti = get_anti(STUB_TR, 'NaN')
    
    else:
        distances_plus = [(tr, tr.start-center) for tr in tr_plus]
        distances_minus = [(tr, center-tr.stop+1) for tr in tr_minus]
       
        distances_plus = [x for x in distances_plus if x[1]>-1*inside]
        distances_minus = [x for x in distances_minus if x[1]>-1*inside]
        
        if(distances_plus):
            transcript_plus, mindistance_plus = min(distances_plus, key = lambda x: abs(x[1]))
        else:
            transcript_plus, mindistance_plus = STUB_TR, float("inf")
            
        if(distances_minus):
            transcript_minus, mindistance_minus = min(distances_minus, key = lambda x: abs(x[1]))
        else:
            transcript_minus, mindistance_minus = STUB_TR, float("inf")
            
            
        if(mindistance_plus < mindistance_minus):
            transcript, mindistance = transcript_plus, mindistance_plus
            anti = get_anti(transcript_minus, mindistance_minus)
        else:
            transcript, mindistance = transcript_minus, mindistance_minus
            anti = get_anti(transcript_plus, mindistance_plus)
        

        if(abs(mindistance) <= maxd):
            gtype = 'upstream'
        else:
            pairs = [(tr, center-tr.start, tr.stop - center -1) for tr in tr_plus + tr_minus]
            pairs = [x for x in pairs if x[1]>=0 and x[2]>=0];
            if(pairs):
                gtype = 'gene'
                transcript = pairs[0][0]
                if(transcript.strand == '+'):
                    mindistance = -1*pairs[0][1]
                else:
                    mindistance = -1*pairs[0][2]
            else:
                gtype = 'intergenic'
                                                  
                
    atg = mindistance + float(transcript.attrs['distance'])
    if(str(atg) != "nan"):
        atg = "%d" % atg
    
    attrs = [("Name", interval.name), ("annotation", transcript.attrs['annotation']), ("function", transcript.attrs['function']), ("gene", transcript.name), ("genesymbol", transcript.attrs['genesymbol']), ("cg", transcript.attrs.get('cg', 'unknown')), ("tss", mindistance), ("atg", atg), ("gtype", gtype)] + anti
    
    if(if_bed):
        return construct_gff_interval( interval.chrom, interval.start, interval.stop, 'annotated', score=interval.score, strand=transcript.strand, source='annotate.py', frame='.', attrs=attrs )
    else:
        for attr_name, attr_value in attrs:
            interval.attrs[attr_name] = str(attr_value);
        interval.strand = transcript.strand
        return interval;


######################################################################################################
### Run analyses ###

comments = read_comments(args.path);
intervals = BedTool(args.path);
if_bed = args.path.split(".")[-1] == 'bed'
tr_dict = defaultdict(lambda: ([], []))

for tr in BedTool(args.transcripts):
    if(tr.strand == '+'):
        tr_dict[tr.chrom][0].append(tr)
    else:
        tr_dict[tr.chrom][1].append(tr) 

annintervals = [annotate_position(interval, tr_dict, args.maxd, args.inside) for interval in intervals]

        
######################################################################################################
### Output Results ###
for comment in comments:
    print(comment)
for interval in annintervals:
    sys.stdout.write(str(interval));
    
  
###################################################################################################################
### DRAWING SECTION ###

def fix_density(ax, bins, scale):
    factor = bins[1]-bins[0]
    ypos = [float(x) for x in ax.get_yticks().tolist()]
    yvals = [float(x)*factor for x in ax.get_yticks().tolist()]
    
    step = scale/yvals[1]*ypos[1]
    newvals = [x for x in np.arange(0, max(yvals), scale)]
    newpos = np.arange(0, step*len(newvals), step)

    #print(sum(ylabels))
    ax.set_yticks(newpos)
    ax.set_yticklabels(["%d" % (x*100) for x in newvals])
    

### Plot a piechart of target types ###
fontsize = 16
colors = ['gold', 'lightblue', 'gray'];
labels = ['upstream', 'gene', 'intergenic']
sizes = Counter([x.attrs['gtype'] for x in annintervals])
sizes = [sizes[x] for x in labels]
labels = ["%s (%d)" % x for x in zip(labels, sizes)]

plt.figure(figsize=(12,9))
plt.pie(sizes, explode=None, labels=labels, colors=colors, shadow=False, startangle=30, pctdistance=0.8, autopct="%1.1f%%", textprops={'fontsize': fontsize})
plt.axis('equal')

plt.title("Upstream: Peak is up to %d bp upstream and %d downstream to the closest TSS" % (args.maxd, args.inside), fontsize = fontsize)
plt.savefig(os.path.join(args.outdir, "interval_genomic_types.%s") % args.format, format = args.format)


### Plot a distribution of TSS distances ###

scores = [float(x.attrs['tss']) for x in annintervals if x.attrs['tss'] != 'nan']
scores.sort();
selected_scores = [x for x in scores if x <= 300 and x >= -100]
#print(min(scores))

fig, axes = plt.subplots(ncols=2, figsize = (22, 7), frameon=False)
fig.tight_layout(rect=[0.05, 0.1, 1, 1])
fig.subplots_adjust(wspace = 0.2)
for data, ax in zip([scores, selected_scores], axes):
    _, bins, _ = ax.hist(data, bins = 20, density = True)
    
    ax.set_xlabel('TSS distance', fontsize=fontsize)
    ax.set_ylabel('Fraction [%]', fontsize=fontsize)    
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fix_density(ax, bins, 0.05)
plt.savefig(os.path.join(args.outdir, "tss_hist.%s") % args.format, format = args.format)



    
    
    






























