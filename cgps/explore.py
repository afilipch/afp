#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Explores various aspects of cgps binding'''

import argparse
import os
import sys
#import scipy
import numpy as np;
import pandas as pd;
from pybedtools import BedTool
from Bio import SeqIO

#from afbio.filters import dsk, usk, trk;
#from afbio.peaks import convolute


parser = argparse.ArgumentParser(description='Explores various aspects of cgps binding');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the binding peaks");
parser.add_argument('--coverage', nargs = '?', required=True, type = str, help = "Path to the binding coverage track");
parser.add_argument('--gctrack', nargs = '?', required=True, type = str, help = "Path to the GC track");
parser.add_argument('--gcdrops', nargs = '?', required=True, type = str, help = "Path to the GC drops");
parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format");


parser.add_argument('--plot', nargs = '?', type = str, help = "Path for the output coverage plot");

args = parser.parse_args();
