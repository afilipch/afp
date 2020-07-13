#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Adds motif presense to the all-in table'''


import argparse
import sys
import os

import pandas as pd
import numpy as np

from sklearn.decomposition import NMF
from sklearn.preprocessing import MinMaxScaler



parser = argparse.ArgumentParser(description='Adds motif presense to the all-in table')#, formatter_class = argparse.RawTextHelpFormatter);
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the master table");
parser.add_argument('--processing', nargs = '?', default=False, const=True, type = bool, help = "If set, new correct table is created for further use")
#parser.add_argument('--fimo', nargs = '?', required=True, type = str, help = "Path to the fimo output, tsv file")
args = parser.parse_args();


def normalize_to_sum(l):
    #print(l)
    norma = sum(l)
    if(norma):
        return [x/norma for x in l]
    else:
        print(l)
        return l
    
    
def normalize_with_zscore(l):
    mean = np.mean(l)
    std = np.std(l);
    zscores = [ (x-mean)/std for x in l ]
    print(min(zscores), max(zscores));
    
def normalize_with_zscore_nonzero(l):
    nzero = [x for x in l if x]
    normalize_with_zscore(nzero)



if(args.processing):
    df = pd.read_csv(args.path, index_col = 1)
    del(df['Set ID'])
    df = df.fillna(0)
    df = df.apply(lambda x: x.str.replace(',', '.').replace('#DIV/0!', 0).astype('float'), axis=1)
    df = df.fillna(0)
    sys.stdout.write(df.to_csv(sep='\t'))
    sys.exit()
    
    
df = pd.read_csv(args.path, index_col = 0, sep='\t')
#Remove all-zero rows
df = df.loc[(df!=0).any(axis=1)]
data = df.to_numpy()
cond_names = list(df.index.values)
gene_names = list(df.columns.values)
 
#Normalize experiments
data = np.apply_along_axis( normalize_to_sum, axis=1, arr=data )
np.apply_along_axis( normalize_with_zscore_nonzero, axis=0, arr=data )

#print(data)
sys.exit()
#Normalize gene expression
scaler = MinMaxScaler()
data = scaler.fit_transform(data);
#print(len(data.max(axis=0)))


sys.exit()
#fit with NMF
model = NMF(n_components=100, verbose=False)
W = model.fit_transform(data)
print(model.reconstruction_err_)
#print(W.shape)
H = model.components_
#for l in H:
    #print(len(l))






