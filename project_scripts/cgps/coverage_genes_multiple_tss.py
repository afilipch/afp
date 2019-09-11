'''Draws coverage plots for the genes with multiple tss'''
import argparse
import os
import sys
import numpy as np;
from collections import defaultdict
from pybedtools import BedTool, Interval
from Bio import SeqIO
import matplotlib.pyplot as plt;
from matplotlib.patches import Rectangle
import copy

from afbio.sequencetools import get_at_content, sliding_window, coverage2dict




parser = argparse.ArgumentParser(description='Explores positions of binding peaks on a genome');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks, gff format");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts coordinates, gff file");
parser.add_argument('--cds', nargs = '?', required=True, type = str, help = "Path to the cds coordinates, gff/bed file");
parser.add_argument('--coverage', nargs = '?', required=True, type = str, help = "Path to the normalized coverage track, bed format");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--length', nargs = '?', default=200, type = int, help = "Length of the upstream/downstream regions");
parser.add_argument('--shorter', nargs = '?', default=0, type = int, help = "If set, only this number of nucleotides in upstream CDS will be shown");
parser.add_argument('--zscore', nargs = '?', default=2, type = int, help = "Z-score threshold for the binding peaks");
parser.add_argument('--at_length', nargs = '?', default=30, type = int, help = "Length of an AT-rich motif");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory")
args = parser.parse_args();


def get_gene_at_profile(interval, genome, window):
    flank = window//2
    seq = genome[interval.chrom][interval.start-flank:interval.stop+flank].seq
    seq = str(seq.upper());
    return [get_at_content(x) for x in sliding_window(seq, window+1)]
    



def compile_transcripts(transcripts, cds_dict, flank_length, shorter):
    genedict = defaultdict(list);
    for t in transcripts:
        genedict[t.name].append(t);
    
    for name, gene_transcripts in genedict.items():
        if(len(gene_transcripts)>1 and len(gene_transcripts)<5):
            gene = Interval(gene_transcripts[0].chrom, min([x.start for x in gene_transcripts]) - flank_length, max([x.stop for x in gene_transcripts]) + flank_length, name, '0', gene_transcripts[0].strand);
            cds = cds_dict[name]
            blocks = [];
            if(gene.strand == '+'):
                curstop = cds.start
                for t in sorted(gene_transcripts, key=lambda x: x.start, reverse = True):
                    blocks.append((t.start, curstop))
                    curstop = t.start;
            else:
                curstart = cds.stop
                for t in sorted(gene_transcripts, key=lambda x: x.stop):
                    blocks.append((curstart, t.stop))
                    curstart = t.stop;
                    
            if(shorter):
                if(gene.strand == '+'):
                    gene = Interval(gene.chrom, gene.start, cds.start+shorter, gene.name, '0', gene.strand)
                    cds = Interval(cds.chrom, cds.start, min(cds.start+shorter, cds.stop), cds.name, '0', cds.strand)
                else:
                    gene = Interval(gene.chrom, cds.stop - shorter, gene.stop, gene.name, '0', gene.strand)
                    cds = Interval(cds.chrom, max(cds.start, cds.stop-shorter), cds.stop, cds.name, '0', cds.strand)  
                
            
            yield (gene, cds, blocks)
            
            
            
            #for t, b in zip(gene_transcripts, blocks):
                #sys.stdout.write(str(t))
                #print(b)
            #print()
            #sys.stdout.write(str(cds))
            #sys.stdout.write(str(gene))
            #print()
            #print("-"*200);
            


            
coverage = coverage2dict(args.coverage)
genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
transcripts = BedTool(args.transcripts);
cds_dict = dict([(x.name, x) for x in BedTool(args.cds)]);
phages = BedTool(args.phages);
regions = BedTool(args.path)
regions = BedTool([x for x in regions if float(x.score)>args.zscore])
OFFSET = regions.field_count()

phaged = transcripts.intersect(b=phages, u =True, f = 0.5)
selected_genes = set(); 
compiled_transcript_list = [x for x in compile_transcripts(phaged, cds_dict, args.length, args.shorter)]
print(len(compiled_transcript_list))
genes = BedTool([x[0] for x in compiled_transcript_list])
for gene in genes.intersect(b=regions, F=0.5, u=True):
    selected_genes.add(gene.name);
compiled_transcript_list = [x for x in compiled_transcript_list if x[0].name in selected_genes]
print(len(compiled_transcript_list))


def draw(compiled_transcript, coverage, genome, at_length, shorter = False, fontsize=28, linewidth=5):
    gene, cds, blocks = compiled_transcript;          
    
    
    local_coverage = coverage[gene.chrom][gene.start:gene.end]
    utr_colors = ['darkblue', 'lightblue', 'blue', 'cyan']
    x_range = range(len(gene))
    at_profile = get_gene_at_profile(gene, genome, args.at_length)
    

    fig, ax1 = plt.subplots(figsize=(20,8))
    ax2 = ax1.twinx()
    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])

    ax1.set_xlabel("gene'[nt]", fontsize=fontsize)
    ax1.set_ylabel('CgpS Coverage [nomalized]', fontsize=fontsize)
    ax2.set_ylabel('AT content', fontsize=fontsize)
    ax1.tick_params(axis='both', labelsize=fontsize, top=False, right=False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    for axis in ['bottom','left','right']:
        ax1.spines[axis].set_linewidth(linewidth)
        
    bottom = max(local_coverage)*(-0.1)     
    ax1.set_ylim(bottom=bottom, top = max(local_coverage)*1.1)
    ax2.set_ylim(bottom=min(at_profile)-0.1, top = max(at_profile)*1.1)
    ax2.plot(x_range, at_profile, color='purple', linewidth=linewidth/2, linestyle='dashed');
    ax1.plot(x_range, local_coverage, color='orange', linewidth=linewidth);
    
    
    #Draw rectangles
    height = -0.8*bottom
    ax1.add_patch(Rectangle((cds.start-gene.start, bottom), len(cds), height, color='green'))
    for block, color in zip(blocks, utr_colors):
        ax1.add_patch(Rectangle((block[0]-gene.start, bottom), block[1]-block[0], height, color=color))
        
    
    #plt.show()
    
    #for profile, color, label in zip(transcript_profiles, colors, labels):
        #ax.plot(x_range, profile, color = color, linewidth=linewidth, label=label)

    #fig.legend(loc=(0.15, 0.86), frameon=False, fontsize=fontsize, ncol = 2)
    if(shorter):
        plt.savefig(os.path.join(args.outdir, "%s_upstream_profile.%s"  %  (gene.name,args.format) ) , format = args.format)
    else:
        plt.savefig(os.path.join(args.outdir, "%s_profile.%s"  %  (gene.name,args.format) ) , format = args.format)
    plt.clf()
    plt.close()

for ct in compiled_transcript_list:
    draw(ct, coverage, genome, args.at_length, shorter = bool(args.shorter), fontsize=28, linewidth=5);
























