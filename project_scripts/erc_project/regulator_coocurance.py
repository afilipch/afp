'''Checks co-occurrence of regulators in among multiple bacterial genomes'''

import argparse
import os
import sys
from itertools import combinations
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt









parser = argparse.ArgumentParser(description='Checks co-occurrence of regulators in among multiple bacterial genomes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the occurrence table, csv format");
#parser.add_argument('--maxlength', nargs = '?', default=10e10, type = int, help = "Maximal length allowed for prophages");
args = parser.parse_args();

def normalized_dot(l1, l2):
    norma = sum(l1)*sum(l2)/len(l1)
    return np.dot(l1, l2)/norma

df = pd.read_csv('/hdd/projects/erc_project/data/filtered_regulators.csv', index_col = 0, delimiter = ',')
data = df.to_numpy()
genome_names = list(df.index.values)
regulator_names = list(df.columns.values)

occurrence_list = []
for i in range(data.shape[1]):
    occurrence_list.append([ 1 if x else 0 for x in data[:,i] ])

signal = []
for i, j in combinations(range(len(occurrence_list)), 2):
    signal.append((i, j, normalized_dot(occurrence_list[i], occurrence_list[j])))
                               
noise = []
for _ in range(1000):
    for i, j in combinations(range(data.shape[1]), 2):
        a = occurrence_list[i]
        b = occurrence_list[j]
        noise.append(normalized_dot(random.sample(a, len(a)), random.sample(b, len(b))))

p_levels = [99.99, 99.9, 99, 95]
p_vals = [np.percentile(noise, x) for x in p_levels]


for i, j, val in signal:
    for p, level in zip(p_vals, p_levels):
        if(val>p):
            print("\t".join((regulator_names[i], regulator_names[j], "%1.1f" % val, "%1.4f" % (1 - level/100) )))
            break;
        




