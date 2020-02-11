#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Annotates the discovered peaks'''

import argparse
import sys
import os
from collections import defaultdict, Counter
from afbio.sequencetools import coverage2dict


import numpy as np;
import matplotlib.pyplot as plt;
from pybedtools import BedTool

from afbio.pybedtools_af import construct_gff_interval, read_comments


parser = argparse.ArgumentParser(description='Annotates the discovered peaks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genomic regions, gff format");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts regions, gff format");
parser.add_argument('--coverage', nargs = '?', type = str, help = "Path to the coverage track, bed format")
parser.add_argument('--inside', nargs = '?', default=200, type = int, help = "Maximum allowed distance to TSS while inside a gene");
parser.add_argument('--maxd', nargs = '?', default=800, type = int, help = "Maximum allowed distance to TSS");
parser.add_argument('--covlimit', nargs = '?', default=1, type = float, help = "Peak boundaries are set when coverage drops to --covlimit")
parser.add_argument('--foldlimit', nargs = '?', default=4, type = float, help = "Peak boundaries are set when coverage is --foldlimit times less than top coverage")

parser.add_argument('--format', nargs = '?', default='svg', type = str, help = "Plot format, svg by default");
parser.add_argument('--outdir', nargs = '?', required=False, type = str, help = "Path to the output plot directory")
args = parser.parse_args();


if(os.stat(args.path).st_size == 0):
    sys.exit("###annotate\nInput file is empty, empty output is produced\n")



STUB_TR = construct_gff_interval( "unknown", 0, 2, 'fake', score='0', strand="+", source='ff', frame='.', attrs= [("Name", "fake"), ("annotation", "None"), ("function", "None"), ("genesymbol", "fake"), ("distance", "NaN")])


### Annotate genomically ###

def annotate_position(peak, tr_dict, maxd, inside):
    center = int(peak.name)
    if(peak.chrom not in tr_dict):
        mindistance = float("NaN")
        gtype = "intergenic"
        transcript = STUB_TR
    
    else:
        center = int(peak.name)
        tr_plus, tr_minus = tr_dict[peak.chrom]
        distances = [(tr, tr.start-center) for tr in tr_plus]
        distances.extend([(tr, center-tr.stop+1) for tr in tr_minus]);
        distances = [x for x in distances if x[1]>-1*inside]
        transcript, mindistance = min(distances, key = lambda x: abs(x[1]))

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
    
       
    #attrs = [("Name", peak.name), ("annotation", transcript.attrs['annotation']), ("function", transcript.attrs['function']), ("gene", transcript.name), ("genesymbol", transcript.attrs['genesymbol']), ("tss", mindistance), ("atg", atg), ("gtype", gtype)]
    attrs = [("Name", peak.name), ("annotation", transcript.attrs['annotation']), ("function", transcript.attrs['function']), ("gene", transcript.name), ("genesymbol", transcript.attrs['genesymbol']), ("cg", transcript.attrs.get('cg', 'unknown')), ("tss", mindistance), ("atg", atg), ("gtype", gtype)]
    
    if(if_bed):
        return construct_gff_interval( peak.chrom, peak.start, peak.stop, 'annotated', score=peak.score, strand=transcript.strand, source='annotate.py', frame='.', attrs=attrs )
    else:
        for attr_name, attr_value in attrs:
            peak.attrs[attr_name] = str(attr_value);
        peak.strand = transcript.strand
        return peak;



tr_dict = defaultdict(lambda: ([], []))

for tr in BedTool(args.transcripts):
    if(tr.strand == '+'):
        tr_dict[tr.chrom][0].append(tr)
    else:
        tr_dict[tr.chrom][1].append(tr)


comments = read_comments(args.path);
peaks = BedTool(args.path);
if_bed = args.path.split(".")[-1] == 'bed' 


annpeaks = [annotate_position(peak, tr_dict, args.maxd, args.inside) for peak in peaks]

        

### Add Peak Intensities ###

def calculate_area_coverage(peak, coverage, covlimit, foldlimit):
    top = int(peak.name);
    topcov = coverage[top]

    for f, v in enumerate(coverage[top:]):
        if(v<covlimit or v*foldlimit<topcov):
            break;
    for b, v in enumerate(coverage[top-1::-1]):
        if(v<covlimit or v*foldlimit<topcov):
            break;
    return sum(coverage[top-b:top+f])
        
if(args.coverage):
    cov_dict = coverage2dict(args.coverage)
    #cov_list = [coverage2dict(args.coverage, cpos=x) for x in range(3, 7)]
    for peak in annpeaks:
        peak.attrs['area_coverage'] = "%1.1f" % calculate_area_coverage(peak, cov_dict[peak.chrom], args.covlimit, args.foldlimit)
        peak.attrs['topcoverage'] =  "%1.3f" % cov_dict[peak.chrom][int(peak.name)]
        #peak.attrs['other_coverage'] =  ",".join(["%1.3f" % x[peak.chrom][int(peak.name)] for x in cov_list])
        




######################################################################################################
### Output Results ###
for comment in comments:
    print(comment)
for peak in annpeaks:
    sys.stdout.write(str(peak));
    
  
  
######################################################################################################
### Draw Plots ###

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
if(not args.outdir):
    sys.exit()

fontsize = 16
colors = ['gold', 'lightblue', 'gray'];
labels = ['upstream', 'gene', 'intergenic']
sizes = Counter([x.attrs['gtype'] for x in annpeaks])
sizes = [sizes[x] for x in labels]
labels = ["%s (%d)" % x for x in zip(labels, sizes)]

plt.figure(figsize=(12,9))
plt.pie(sizes, explode=None, labels=labels, colors=colors, shadow=False, startangle=30, pctdistance=0.8, autopct="%1.1f%%", textprops={'fontsize': fontsize})
plt.axis('equal')

plt.title("Upstream: Peak is up to %d bp upstream and %d downstream to the closest TSS" % (args.maxd, args.inside), fontsize = fontsize)
plt.savefig(os.path.join(args.outdir, "peak_genomic_types.%s") % args.format, format = args.format)
plt.close()


### Plot a scatter of peaks topcoverage/area_coverage ###
if(args.coverage and len(annpeaks)>20):
    fontsize=24
    linewidth = 5
    xvals = np.array([float(x.attrs['topcoverage']) for x in annpeaks])
    yvals = np.array([float(x.attrs['area_coverage']) for x in annpeaks])

    fig, ax = plt.subplots(figsize=(16,9))
    #plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])

    ax.set_xlabel('Top coverage [avg]', fontsize=fontsize)
    ax.set_ylabel('Area coverage [avg]', fontsize=fontsize)    
    ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.scatter(xvals, yvals)
    plt.title("Area boundaries: drop to %1.1f in coverage averages or %1.1f folds to top coverage" % (args.covlimit, args.foldlimit), fontsize = fontsize)
    plt.savefig(os.path.join(args.outdir, "peak_coverage_scores.%s") % args.format, format = args.format)
    plt.close()



### Plot a distribution of peaks topcoverage ###

    scores = [float(x.attrs['topcoverage']) for x in annpeaks]
    scores.sort();
    selected_scores = scores[:int(len(scores)*0.75)]
    #print(min(scores))
    
    fig, axes = plt.subplots(ncols=2, figsize = (22, 7), frameon=False)
    fig.tight_layout(rect=[0.05, 0.1, 1, 1])
    fig.subplots_adjust(wspace = 0.2)
    for data, ax in zip([scores, selected_scores], axes):
        _, bins, _  = ax.hist(data, bins = 20, range = (0, max(data)), density = True)
        ax.set_xlabel('Top coverage', fontsize=fontsize)
        ax.set_ylabel('Fraction [%]', fontsize=fontsize)    
        ax.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        fix_density(ax, bins, 0.05)
    
    plt.savefig(os.path.join(args.outdir, "topcoverage_hist.%s") % args.format, format = args.format)
    plt.close()
    
    
    
### Plot a distribution of peaks TSS distances ###

    scores = [float(x.attrs['tss']) for x in annpeaks if x.attrs['tss'] != 'nan']
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
    plt.close()





























