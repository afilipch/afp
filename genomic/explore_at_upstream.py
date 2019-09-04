#! /usr/bin/python
'''Analyses AT content relative to the transcripts and phages positions on a genome'''

import argparse
import sys
import os
from collections import defaultdict
import matplotlib.pyplot as plt;
import numpy as np;

from Bio import SeqIO;
from pybedtools import BedTool;

from afbio.sequencetools import get_at_content, sliding_window, transcript_upstream
from afbio.numerictools import CDF, lists2thresholds




parser = argparse.ArgumentParser(description='Analyses AT content relative to the transcripts and phages positions on a genome');
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta file");
parser.add_argument('--transcripts', nargs = '?', required=True, type = str, help = "Path to the transcripts coordinates, bed/gff file");
parser.add_argument('--phages', nargs = '?', required=True, type = str, help = "Path to the phages coordinates, bed file");
parser.add_argument('--length', nargs = 2, default=(10, 50), type = int, help = "Length range of the segments upstream");
parser.add_argument('--exact_distance', nargs = '?', default= 70, type = int, help = "Distance (to the transcript start) of the segments upstream for the exact search");
parser.add_argument('--distances', nargs = 2, default=(0, 70), type = int, help = "Distances (min and max distance) of the segments upstream for the max search");
parser.add_argument('--format', nargs = '?', default='png', type = str, help = "Plot format, png by default");
parser.add_argument('--outdir', nargs = '?', required=True, type = str, help = "Path to the output directory");

#parser.add_argument('--tss', nargs = '?', required=True, type = str, help = "Path to the TSS annotation file");
#parser.add_argument('--stranded', nargs = '?', const=True, default=False, type = bool, help = "If set, the data considered to be stranded");
args = parser.parse_args();


def transcript2upstream(interval, genome,  distance, length):
    
    seq = None
    if(interval.strand == '+'):
        start = interval.start - length - distance;
        end = interval.start;
        if(start >= 0 ):
            seq = str(genome[interval.chrom][start:end].seq.reverse_complement().upper())
    elif(interval.strand == '-'):
        start = interval.stop
        end = interval.stop + length + distance;
        chrom = genome[interval.chrom]
        if(end<len(chrom)):
            seq = str(chrom[start:end].seq.upper())
    return seq


def exact_at_mindistance_length(sequences, maxdistance, length_range):
    #Distance is row, length is column
    arr = np.zeros((maxdistance, len(list(length_range))))
    for i, distance in enumerate(range(maxdistance)):
        for j, length in enumerate(length_range):
            arr[i, j] = np.mean([get_at_content(x[distance: distance+length]) for x in sequences])
    return arr;


def max_at_maxdistance_length(sequences, max_distance, length_range):
    
    tail = max(length_range) + 1
    d_range = range(tail, max_distance);
    arr = np.zeros((   len(list(d_range)), len(list(length_range))   ))
                   
    for i, d in enumerate(d_range):
        local_sequences = [x[:d] for x in sequences]
        for j, length in enumerate(length_range):
            temp_list =  []
            for seq in local_sequences:
                temp_list.append(max([get_at_content(x) for x in sliding_window(seq, length)]))
            arr[i,j] = np.mean(temp_list)
    return arr


def max_at_mindistance_length(sequences, max_distance, length_range):
    
    tail = max(length_range)
    d_range = range(0, max_distance-tail);
    arr = np.zeros((   len(list(d_range)), len(list(length_range))   ))
                   
    for i, d in enumerate(d_range):
        local_sequences = [x[d:] for x in sequences]
        for j, length in enumerate(length_range):
            temp_list =  []
            for seq in local_sequences:
                temp_list.append(max([get_at_content(x) for x in sliding_window(seq, length)]))
            arr[i,j] = np.mean(temp_list)
    return arr




#def length_dependence(sequences, min_distance, max_distance, length_range):
    #local_sequences = [x[min_distance:max_distance] for x in sequences]
    #arr = []
    #for length in length_range:
        #temp_list =  []
        #for seq in local_sequences:
            #temp_list.append(max([get_at_content(x) for x in sliding_window(seq, length)]))
        #arr.append(np.mean(temp_list))
    #return np.array(arr)


#def get_3d_tables(sequences, min_distance, max_distance, length_range):
    #tail = (max_distance-min_distance-max(length_range)-1)//2
    #print(tail)
    #d1_range = range(min_distance, min_distance+tail);
    #d2_range = range(max_distance-tail, max_distance);
    #arr = np.zeros((   len(list(d1_range)), len(list(d2_range)), len(list(length_range))   ))
                   
    #for i, d1 in enumerate(d1_range):
        #for j, d2 in enumerate(d2_range):
            #local_sequences = [x[d1:d2] for x in sequences]
            #for k, length in enumerate(length_range):
                #temp_list =  []
                #for seq in local_sequences:
                    #temp_list.append(max([get_at_content(x) for x in sliding_window(seq, length)]))
                #arr[i,j,k] = np.mean(temp_list)
    #return arr




def _local_score(n1, n2, total):
    spec = (n1)/(n1+n2)
    sens = n1/total
    return spec*sens

def separation_score(l1, l2):
    total = len(l1);
    thr_dict = lists2thresholds(l1, l2, True);
    #print(len([x for x in l2 if x >=800]))
    res = [ (x[0], _local_score(x[1][0], x[1][1], total), x[1][0], x[1][1]) for x in thr_dict.items() ]
    return max(res, key=lambda x: x[1]);
    return sum([x[1] for x in res])
    #for r in res:
        #print(r)
    
    



def get_2d_best_separation(sequences1, sequences2, max_distance, length_range):
    tail = max(length_range)
    d_range = range(tail, max_distance);
    arr = np.zeros((   len(list(d_range)), len(list(length_range))   ))
    data = {}
                   
    for i, d in enumerate(d_range):
        local_sequences1 = [x[:d] for x in sequences1]
        local_sequences2 = [x[:d] for x in sequences2]
        for j, length in enumerate(length_range):
            temp_list1, temp_list2 =  [], []
            for seq in local_sequences1:
                temp_list1.append(int(max( [get_at_content(x) for x in sliding_window(seq, length)] ) * 1000) )
            for seq in local_sequences2:
                temp_list2.append(int(max( [get_at_content(x) for x in sliding_window(seq, length)] ) * 1000) )
                
            threshold, score, passed1, passed2 = separation_score(temp_list1, temp_list2);
            arr[i,j] = score
            data[(i, j)] = (d, length, threshold, passed1, passed2);
            
                
            #arr[i,j] = np.mean(temp_list)
    max_indices = np.unravel_index(np.argmax(arr, axis=None), arr.shape)
    print(data[max_indices])
            
            
            
    
    
    


length_range = range(*args.length)
genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
transcripts = BedTool(args.transcripts);
phages = BedTool(args.phages);


phaged_transcripts = [transcripts.intersect(b=phages, u =True, f = 0.5), transcripts.intersect(b=phages, v =True, f = 0.5)]
list_exact_mindistance_length, list_max_maxdistance_length, list_max_mindistance_length = [], [], []
seq_list = [];
for ptr in phaged_transcripts:
    sequences = [transcript2upstream(x, genome, args.distances[1], args.length[1]) for x in ptr]
    sequences = [x for x in sequences if x]
    list_exact_mindistance_length.append(exact_at_mindistance_length(sequences, args.exact_distance, length_range))
    list_max_maxdistance_length.append(max_at_maxdistance_length(sequences, args.distances[1], length_range))
    list_max_mindistance_length.append(max_at_mindistance_length(sequences, args.distances[1], length_range))
    

    
    
    
 
def get_labels_ticks(arr, d_starts, step):
    res = []
    for d, shape in zip(d_starts, arr.shape):
        q, r = divmod(d, step);
        if(r):
            q += 1;
            r = step -r
        ticks = range(r, shape+1, step);
        labels =[(q+x)*step for x in range(len(ticks))]
        res.append((ticks, labels));
    return res;
 
    
def plot_heatmap(cmatrix, d_starts, name, arrtype, tick_step=5, fontsize=12):
    cmap="RdPu"
    
    aspect = [x/sum(cmatrix.shape) for x in cmatrix.shape]
    aspect = [x*20 for x in aspect[::-1]]
    aspect[0] = aspect[0] + 5;
    
    fig, ax = plt.subplots(figsize=aspect)
    plt.tight_layout(rect=[0, 0, 0.8, 0.9])
    im = ax.imshow(cmatrix, cmap=cmap)
    cbar = ax.figure.colorbar(im, ax=ax, cmap=cmap)
    cbar.ax.set_ylabel("AT content %s" % arrtype, rotation=-90, va="bottom", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    
    (yticks, ylabels), (xticks, xlabels) = get_labels_ticks(cmatrix, d_starts, tick_step)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(xlabels, fontsize=fontsize)
    ax.set_yticklabels(ylabels, fontsize=fontsize)
    ax.set_ylabel("distance", fontsize=fontsize)
    ax.set_xlabel("length", fontsize=fontsize)
    ax.xaxis.set_label_position('top') 
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    ax.tick_params(which="minor", bottom=False, left=False)
    plt.title("%s: %s" % (arrtype, name), fontsize=fontsize*1.3)

    plt.savefig(os.path.join(args.outdir, "%s_%s_at.%s"  %  (name, arrtype, args.format)) , format = args.format)
    plt.clf()
    
    

absolute_list = [list_exact_mindistance_length, list_max_maxdistance_length, list_max_mindistance_length]
names = ["exact_mindistance_length", "max_maxdistance_length", "max_mindistance_length"]
adjustments = [(0, args.length[0]), (args.length[1], args.length[0]), (0, args.length[0])]
arrtypes = ['phage', 'non-phage', 'difference']

for arr_list, name, adj in zip(absolute_list, names, adjustments):
    
    arr_list.append(arr_list[0] - arr_list[1]);
    max_indices = [np.unravel_index(np.argmax(x, axis=None), x.shape) for x in arr_list];
    maxes = [x[0][x[1]] for x in zip(arr_list, max_indices)]
    
    print("Best parameters for the scenario: %s" % name.upper())
    for max_index, maxval in zip(max_indices, maxes):
        print("%d\t%d\t%1.2f" % (max_index[0] + adj[0], max_index[1] + adj[1], maxval))
    print()
    
    for cmatrix, arrtype in zip(arr_list, arrtypes):
        plot_heatmap(cmatrix, adj, name, arrtype, tick_step=5, fontsize=20)
        

#get_labels_ticks(list_exact_mindistance_length[2], (0, args.length[0]), 5)    
#plot_heatmap(list_exact_mindistance_length[2], "AT content difference", (0, args.length[0]), tick_step=5, fontsize=16)
    
    
    
    
#get_2d_best_separation(seq_list[0], seq_list[1], args.distance[1], length_range)





    
    


    
    



    
