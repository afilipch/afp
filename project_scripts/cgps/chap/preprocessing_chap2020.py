#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Renames the raw fastq files for better clarity'''


import argparse
import sys
import os
from os import listdir
from os.path import isfile

from afbio.generators import get_only_files
from pathlib import Path
from shutil import copyfile

##Read configuration
#conf = load_config('chipchap')
#bowtie_settings = conf['bowtie'];



parser = argparse.ArgumentParser(description='Renames the fastq files for better clarity')#, formatter_class = argparse.RawTextHelpFormatter);
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with reads");
parser.add_argument('--table', nargs = '?', type = os.path.abspath, required = True, help = "Path to the sample table, tsv format");
parser.add_argument('--outdir', nargs = '?', type = str, required = True, help = "Path to the output directory");
#parser.add_argument('--paired', nargs = '?', default = False, const=True, type = bool, help = "If set, reads are assumed to be paired-end")
#parser.add_argument('--table', nargs = '?', type = os.path.abspath, help = "Path to the table which connects the read file names to the meaningful names");
args = parser.parse_args();


sample2type = {};

with open(args.table) as f:
    for l in f:
        a = l.strip().split("\t");
        time = a[1].replace(" ", "").replace(".", "")
        sample2type[a[2]] = a[0], time, "chap"
        sample2type[a[3]] = a[0], time, "control"


for cond in set([x[0] for x in sample2type.values()]):
    for type_ in ('chap', 'control'):
        path = os.path.join(args.outdir, "%s_%s" % (cond, type_))
        Path(path).mkdir(parents=True, exist_ok=True)
#sys.exit()

for f in get_only_files(args.path):
    if(f.endswith('fastq')):
        name, mate, _ = os.path.basename(f).split(".")
        cond, time, type_ = sample2type[name]
        path = os.path.join(args.outdir, "%s_%s" % (cond, type_), "%s.%s.fastq" % (time, mate))
        copyfile(f, path)
