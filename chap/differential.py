#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Finds consensus regions for the peaks found in different experiments'''

import argparse
import os
import sys
from collections import defaultdict


import numpy as np;
import pandas as pd;
from pybedtools import BedTool, Interval
import matplotlib.pyplot as plt;

from afbio.pybedtools_af import construct_gff_interval


parser = argparse.ArgumentParser(description='Finds consensus regions for the peaks found in different experiments');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "Path to the annotated peaks of the replicates for two different conditions. Conditions must be split by comma, tha is 'a_r1 a_r2 a_r3, b_r1 b_r2'");
#parser.add_argument('--rep_numbers', nargs = '+', required=True, type = str, help = "Numbers of replicates per condition. If we provide 3 replicates for condition '1' and the 4 replicates for condition '2' then we have to set '--rep_numbers 3 4'")
parser.add_argument('--zscore', nargs = '?', default=5.0, type = float, help = "Mininimum z-score required for a peak to be a seed for the consensus region");
parser.add_argument('--maxd', nargs = '?', default=60, type = int, help = "Maximal distance allowed between top positions of the peaks");
#parser.add_argument('--peakflank', nargs = '?', default=5, type = int, help = "Flank's length around a peak to detect maximum coverage");
args = parser.parse_args()

