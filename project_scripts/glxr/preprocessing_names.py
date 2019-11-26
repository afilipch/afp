#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Renames the raw fastq files for better clarity'''


import argparse
import sys
import os
from os import listdir
from os.path import isfile


##Read configuration
#conf = load_config('chipchap')
#bowtie_settings = conf['bowtie'];



parser = argparse.ArgumentParser(description='Renames the raw fastq files for better clarity')#, formatter_class = argparse.RawTextHelpFormatter);
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with reads");
#parser.add_argument('--package', nargs = '?', type = os.path.abspath, required = True, help = "Path to the afp package");
#parser.add_argument('--outdir', nargs = '?', type = str, required = True, help = "Path to the output directory");
#parser.add_argument('--paired', nargs = '?', default = False, const=True, type = bool, help = "If set, reads are assumed to be paired-end")
#parser.add_argument('--table', nargs = '?', type = os.path.abspath, help = "Path to the table which connects the read file names to the meaningful names");
args = parser.parse_args();

onlyfiles = [x for x in listdir(args.path) if isfile(os.path.join(args.path, x))]
for f in onlyfiles:
    a = f.strip("\'").split("_")
    for el in a:
        if(el in ["R1", "R2"]):
            mate = el;
            break;
    num = a[0].split("-")[0];
    if(len(num)>3):
        num = num[1:];
    print(num)
    os.rename(os.path.join(args.path, f), os.path.join(args.path, "%s_%s.fastq" % (num, mate)))
