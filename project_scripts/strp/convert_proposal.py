
import argparse
import os
import sys
from collections import Counter



parser = argparse.ArgumentParser(description='converts');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to txt file");
#parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Path to the table with proteins from human genome atlas");
args = parser.parse_args()


dictionary = {"viruses": 'neurons', "surface": "p-body", "bacterial": "ancestral", "host": "liver",  "bacteria" : 'ribosome', "peptide": "pseudogene", "attachment" : 'aggregation', "virus": "neuron", "microbiome" : 'transcriptome', "protein" : 'contamination'}



with open(args.path) as f:
    res = f.read()
    for from_, to_ in dictionary.items():
        res = res.replace(from_, to_ + '*')
        
sys.stdout.write(res)
    
    #words = Counter(res.replace('\n', '').replace('\t', '').split(' '))
 
#print('"' + '", "'.join([x[0] for x in words.items() if x[1] > 2]) + '"')
      
#for w, c in words.items():
    #if(c>=3):
        #print (w, c)
