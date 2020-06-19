#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Convert proprietary TOPDOM xml file into more convenient tsv'''

import argparse
import os
import sys
from collections import defaultdict

import xml.etree.ElementTree as ET
import csv


parser = argparse.ArgumentParser(description='Convert proprietary TOPDOM xml file into more convenient tsv');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to TOPDOM xml file");
parser.add_argument('--specie', nargs = '?', type = str, help = "Select entries only for the given specie");
parser.add_argument('--table', required = True, nargs = '?', type = str, help = "Script also creates a table which connects entry name with identical sequences");
args = parser.parse_args()

def create_name(protein, region):
    pname, organism = protein.attrib["ID"].split("_")
    return ":".join(( organism, pname, region.attrib['Begin'], region.attrib['End'] ))

tree = ET.parse(args.path)
root = tree.getroot()
seq2name = defaultdict(list)

with open(args.table, 'w') as f:
    for entry in root:
        for protein in entry:
            if(not args.specie or args.specie == protein.attrib["ID"].split("_")[1]):
                for region in protein:
                    for seq in region:
                        seq2name[seq.text].append(create_name(protein, region));
                        
    for seq, names in seq2name.items():
        print(">%s\n%s" % (names[0], seq))
        f.write("\t".join(names) + "\n")
        
    #break;
        
    #break;
    #break
    #for ell in el:
        #print(ell)

#for el in root.findall('TOPDOM'):
    #print(list(el.items()))
