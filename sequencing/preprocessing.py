#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates makefile for the reads preprocessing'''


import argparse
import sys
import os
from os import listdir
from os.path import isfile

from afbio.config.config import load_config;
from afbio.makefiles import dependence, get_header, get_bowtie_call, get_script, get_bowtie_help

##Read configuration
#conf = load_config('chipchap')
#bowtie_settings = conf['bowtie'];



parser = argparse.ArgumentParser(description='Creates makefile for the reads preprocessing')#, formatter_class = argparse.RawTextHelpFormatter);
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the folder with reads");
parser.add_argument('--package', nargs = '?', type = os.path.abspath, required = True, help = "Path to the afp package");
parser.add_argument('--outdir', nargs = '?', type = str, required = True, help = "Path to the output directory");
parser.add_argument('--paired', nargs = '?', default = False, const=True, type = bool, help = "If set, reads are assumed to be paired-end")
parser.add_argument('--table', nargs = '?', type = os.path.abspath, help = "Path to the table which connects the read file names to the meaningful names");
args = parser.parse_args();

def split_into_pairs(path):
    onlyfiles = [os.path.join(path, x) for x in listdir(path) if isfile(os.path.join(path, x))]
    onlyfiles.sort();
    return list(zip(onlyfiles[::2], onlyfiles[1::2]))

def get_names(table):
    names = {}
    with open(table) as f:
        for l in f:
            a = l.strip().split("\t")
            names[a[0]] = a[1];
    return names


########################################################################################################################
## Main function to create one-sample Makefile
def makefile_local(mates, names, outdir):
    final_files = []
    mlist=[];
    
    for pair in mates:
    # Detect peaks
        rname = os.path.basename(pair[0]).split("_")[0]
        name = names.get(rname)
        if(name):
            input_files = pair
            output_files = [os.path.join(outdir, "%s.%d.fastq" % (name, x)) for x in (1,2)]
            script = get_script('collapse_reads.py', seq_package, arguments={'--output': os.path.join(outdir,name)}, inp = input_files)
            mlist.append(dependence(input_files, output_files, script));  
            final_files.extend(output_files)
        else:
            sys.stderr.write("\n%s name is not found in the provided table\n" % rname)
    
    #Get header and cleaner for the makefile
    sys.stderr.write("\n" + " ".join([x for x in final_files if 'control' not in x]) + "\n\n")
    sys.stderr.write(" ".join([x for x in final_files if 'control' in x]) + "\n\n")
    mlist.insert(0, get_header(final_files))
    mlist.append('clean:\n\techo "nothing to clean."\n');

    return "\n\n".join(mlist)


seq_package = os.path.join(args.package, 'sequencing')
mates = split_into_pairs(args.path)
names = get_names(args.table)
print(makefile_local(mates, names, args.outdir))

