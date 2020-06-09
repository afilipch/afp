#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Convert proprietary TOPDOM xml file into more convenient tsv'''

import argparse
import os
import sys

import xml.etree.ElementTree as ET
import csv


parser = argparse.ArgumentParser(description='Convert proprietary TOPDOM xml file into more convenient tsv');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to TOPDOM xml file");
#parser.add_argument('--table', nargs = '?', required=True, type = str, help = "Path to the table with proteins from human genome atlas");
args = parser.parse_args()

tree = ET.parse(args.path)
root = tree.getroot()
for el in root:
    print(el.attrib['Source'])
    #break
    #for ell in el:
        #print(ell)

#for el in root.findall('TOPDOM'):
    #print(list(el.items()))
