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
parser.add_argument('--distance', nargs = 2, default=(0, 70), type = int, help = "Distance (to the transcript start) range of the segments upstream");
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


def get_at_tables(sequences, distance_range, length_range):
    #Distance is row, length is column
    arr = np.zeros((len(list(distance_range)), len(list(length_range))))
    for i, distance in enumerate(distance_range):
        for j, length in enumerate(length_range):
            cur_at = np.mean([get_at_content(x[distance: distance+length]) for x in sequences])
            arr[i, j] = cur_at;
    return arr;

def length_dependence(sequences, min_distance, max_distance, length_range):
    local_sequences = [x[min_distance:max_distance] for x in sequences]
    arr = []
    for length in length_range:
        temp_list =  []
        for seq in local_sequences:
            temp_list.append(max([get_at_content(x) for x in sliding_window(seq, length)]))
        arr.append(np.mean(temp_list))
    return np.array(arr)


def get_3d_tables(sequences, min_distance, max_distance, length_range):
    tail = (max_distance-min_distance-max(length_range)-1)//2
    print(tail)
    d1_range = range(min_distance, min_distance+tail);
    d2_range = range(max_distance-tail, max_distance);
    arr = np.zeros((   len(list(d1_range)), len(list(d2_range)), len(list(length_range))   ))
                   
    for i, d1 in enumerate(d1_range):
        for j, d2 in enumerate(d2_range):
            local_sequences = [x[d1:d2] for x in sequences]
            for k, length in enumerate(length_range):
                temp_list =  []
                for seq in local_sequences:
                    temp_list.append(max([get_at_content(x) for x in sliding_window(seq, length)]))
                arr[i,j,k] = np.mean(temp_list)
    return arr


def get_2d_tables(sequences, max_distance, length_range):
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


def _local_score(n1, n2, total):
    spec = (n1)/(n1+n2)
    sens = n1/total
    return (spec**3)*(sens**2)

def separation_score(l1, l2):
    total = len(l1);
    thr_dict = lists2thresholds(l1, l2, True);
    #print(len([x for x in l2 if x >=800]))
    res = [ (x[0], _local_score(x[1][0], x[1][1], total), x[1][0], x[1][1]) for x in thr_dict.items() ]
    return max(res, key=lambda x: x[1]);
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
            
            
            
    
    
    


distance_range = range(*args.distance)
length_range = range(*args.length)
genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
transcripts = BedTool(args.transcripts);
phages = BedTool(args.phages);


phaged_transcripts = [transcripts.intersect(b=phages, u =True, f = 0.5), transcripts.intersect(b=phages, v =True, f = 0.5)]
seq_list = [];
for ptr in phaged_transcripts:
    sequences = [transcript2upstream(x, genome, args.distance[1], args.length[1]) for x in ptr]
    sequences = [x for x in sequences if x]
    seq_list.append(sequences)
get_2d_best_separation(seq_list[0], seq_list[1], args.distance[1], length_range)



#at_tables = [];
#for ptr in phaged_transcripts:
    #sequences = [transcript2upstream(x, genome, args.distance[1], args.length[1]) for x in ptr]
    #sequences = [x for x in sequences if x]
    ##at_tables.append(get_at_tables(sequences, distance_range, length_range))
    ##at_tables.append(length_dependence(sequences, args.distance[0], args.distance[1], length_range))
    ##at_tables.append(get_3d_tables(sequences, args.distance[0], args.distance[1], length_range))
    #at_tables.append(get_2d_tables(sequences, args.distance[1], length_range))
    
    


##print(at_tables)
#at_tables.append(at_tables[0] - at_tables[1]);
#max_indices = [np.unravel_index(np.argmax(x, axis=None), x.shape) for x in at_tables];
#print(max_indices)
#maxes = [x[0][x[1]] for x in zip(at_tables, max_indices)]
#print(maxes);


    
    


    
    



    
