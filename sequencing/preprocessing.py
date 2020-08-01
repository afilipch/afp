#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates makefile for the reads preprocessing'''


import argparse
import sys
import os
from collections import defaultdict

from afbio.config.config import load_config;
from afbio.makefiles import dependence, get_header, get_bowtie_call, get_script, get_bowtie_help

from afbio.generators import get_only_files
from time import gmtime, strftime

parser = argparse.ArgumentParser(description='Creates makefile for the reads preprocessing')#, formatter_class = argparse.RawTextHelpFormatter);
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with reads");
parser.add_argument('--package', nargs = '?', type = os.path.abspath, required = True, help = "Path to the afp package");
parser.add_argument('--outdir', nargs = '?', type = str, required = True, help = "Path to the output directory");
parser.add_argument('--paired', nargs = '?', default = False, const=True, type = bool, help = "If set, reads are assumed to be paired-end")
parser.add_argument('--table', nargs = '?', required = True, type = os.path.abspath, help = "Path to the table which connects the read file names to the meaningful names");
parser.add_argument('--log', nargs = '?', required = True, type = os.path.abspath, help = "Path to the log file");
args = parser.parse_args();






########################################################################################################################
## Main function to create one-sample Makefile
def makefile_local(name2sample, outdir, paired):
    final_files = []
    mlist=[];
    
    if(paired):
        for name, pair in name2sample.items():
            input_files = pair
            output_files = [os.path.join(outdir, "%s.%d.fastq" % (name, x)) for x in (1,2)]
            script = get_script('collapse_paired.py', seq_package, arguments={'--output': os.path.join(outdir,name), '--log': args.log}, inp = input_files)
            mlist.append(dependence(input_files, output_files, script));  
            final_files.extend(output_files)

    else:
        for name, pair in name2sample.items():
            input_files = pair
            output_files = os.path.join(outdir, "%s.fastq" % name)
            script = get_script('collapse_single.py', seq_package, arguments={'--output': os.path.join(outdir, name), '--log': args.log}, inp = input_files)
            mlist.append(dependence(input_files, output_files, script));  
            final_files.append(output_files)
           
    
    mlist.insert(0, get_header(final_files))
    mlist.append('clean:\n\techo "nothing to clean."\n');

    return "\n\n".join(mlist)


seq_package = os.path.join(args.package, 'sequencing')

sample2name = {}
with open(args.table) as f:
    for l in f:
        a = l.strip().split("\t")
        sample2name[a[0]] = a[1:];

if(args.paired):
    name2sample = defaultdict(lambda: [None]*2)
    for sample in get_only_files(args.path):
        a = sample2name.get(os.path.basename(sample))
        if(a):
            name2sample[a[0]][int(a[1])-1] = sample
        else:
            sys.stderr.write("Sample %s was not found in the provided table\n" % sample)
else:
    name2sample = {}
    for sample in get_only_files(args.path):
        a = sample2name.get(os.path.basename(sample))
        if(a):
            name2sample[a[0]] = sample
        else:
            sys.stderr.write("Sample %s was not found in the provided table\n" % sample)


print(makefile_local(name2sample, args.outdir, args.paired))
with open(args.log, 'w') as f:
    timestr = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    f.write("Time created: %s\n\nProject call: python %s\n" % (timestr, " ".join(sys.argv)) );

