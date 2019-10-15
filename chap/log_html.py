import sys
import os
import dominate
import argparse
from dominate.tags import *
from dominate.util import raw


from afbio.html.methods import add_ucsc
from math import log

parser = argparse.ArgumentParser(description='Generates html report of CHAP analyses by chipchap.py');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the log file, txt format");
#parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format");
#parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
#parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
#parser.add_argument('--flank', nargs = '?', default=60, type = int, help = "Peak plank length");
args = parser.parse_args();

name = ".".join(os.path.basename(args.path).split(".")[:-1])
print(name)
