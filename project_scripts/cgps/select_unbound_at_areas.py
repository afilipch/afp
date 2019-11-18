'''Very simple selection of the unbound at_rich areas'''
import argparse
import os
import sys
from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Very simple selection of the unbound at_rich areas');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the at rich areas, gff format");
parser.add_argument('--threshold', nargs = '?', required=True, type = float, help = "AT content threshold");
args = parser.parse_args();

for interval in BedTool(args.path):
    if(interval.attrs['type'] == 'out' and float(interval.score)>=args.threshold):
        interval.start -= 25
        interval.stop += 25
        sys.stdout.write(str(interval))
