#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates makefile and directory structure for chip/chap seq analysis'''

from time import gmtime, strftime
import argparse
import sys
import os
from os import listdir
from os.path import isfile, isdir

from afbio.config.config import load_config;
from afbio.makefiles import dependence, get_header, get_bowtie_call, get_script, get_bowtie_help

#Read configuration
conf = load_config('chipchap')
bowtie_settings = conf['bowtie'];
coverage_settings = conf['coverage']
annotation_settings = conf['annotation']
region_settings = conf['region']
control_settings = conf['control']

#Bowtie2 arguments preliminary parsing
bowtie_help_str = "[%s]" % " ".join(["%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in bowtie_settings.items()])


parser = argparse.ArgumentParser(description='Creates makefile and directories for chiflex project')#, formatter_class = argparse.RawTextHelpFormatter);
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = os.path.abspath, help = "Path to the project folder. If folder does not exist, it will be created");
parser.add_argument('--package', nargs = '?', required=True, type = os.path.abspath, help = "Path to the afp package");
parser.add_argument('--reads', nargs = '+', required=True, type = os.path.abspath, help = "Path to sequencing reads. fastq/fasta file. Paired-end reads must be provided consecutively");
parser.add_argument('--index', nargs = '?', required=True, type = os.path.abspath, help = "Path to the mapping reference bowtie2 indices");
parser.add_argument('--genome', nargs = '?', required=True, type = os.path.abspath, help = "Path to the mapping references in fasta format");
parser.add_argument('--control', nargs = '+', type = os.path.abspath, help = "Path to the genomic control fasta files used for the adjustment of multiple dna copies. The order of the files must be the same as for --reads option. If for some --reads files corresponding controls are absent, \'none\' arguments MUST be provided in the corresponding places");

#Options for the mapping result postprocessing
parser.add_argument('--collapsed', nargs = '?', default = False, const=True, type = bool, help = "If set, reads are assumed collapsed with collapse.pl script. Read count appendix of the read id will be used to calculate read support of the interactions")
parser.add_argument('--paired', nargs = '?', default = False, const=True, type = bool, help = "If set, reads are assumed to be paired-end")
#Options for the output control
parser.add_argument('--only_makefile', nargs = '?', default = False, const = True, type = bool, help = "If set, a new makefile is created, but not folders");

#Options controlling postprocessing
parser.add_argument('--multi', nargs = '?', default = False, const = True, type = bool, help = "If set, consensus regions will be generated and explored");
parser.add_argument('--annotation', nargs = '?', default=False, type = os.path.abspath, help = "Path to genbank annotation file in gff format.");
parser.add_argument('--ucsc', nargs = '?', default="stub", type = str, help = "Name of the UCSC session for the experiment");
parser.add_argument('--name', nargs = '?', default="stub", type = str, help = "Name of the project");


#bowtie2  and filtering options
parser.add_argument('--bowtie_args', nargs = '+', default = [], type = str, help = "Bowtie settings for the first round of mapping. For example, if one wants to set \'-p 4\', use \'--local\' alignment mode, but not \'--norc\' option then \'p=4 local=True norc=False\' should be provided. Given attributes replace default(for Chiflex, NOT for Bowtie) ones. Default settings for the modes are: %s" % bowtie_help_str)
parser.add_argument('--reads_format', nargs = '?', default = 'U', choices = ['U', 'f'], type = str, help = "Reads format: 'U' -> fastq, 'f' -> fasta");

parser.add_argument('--filtering', nargs = '?', default="normal", choices=['loose', 'normal', 'strict'], type = str, help = "Strength of filtering, can be: 'loose', 'normal' or 'strict'");
parser.add_argument('--ambiguous', nargs = '?', default=0, const=1, type = int, help = "If, set ambiguous mappings will be also counted");

#Performance option
parser.add_argument('--threads', nargs = '?', default = 1, type = int, help = "Number of threads to use");

args = parser.parse_args();



#######################################################################################################################
# Process input options

peak_filtering_settings = conf['peak_filtering_%s' % args.filtering]

#######################################################################################################################
# Set paths' constants

rnaseq_package = os.path.join(args.package, 'rnaseq')
mapping_package = os.path.join(args.package, 'mapping')
html_lib = os.path.join(args.package, 'afbio', 'html')


#######################################################################################################################
# Create folders

project_path = args.path
folders = ['sam', 'coverage', 'log', 'transcripts', 'makefiles', 'differential', 'ucsc']


while(not args.only_makefile):
    try:
        os.makedirs(project_path);
        for folder in folders:
            os.mkdir(os.path.join(project_path, folder));
        break;
    except:
        answer = input("\nProject directory \'%s\' is currently exists, please type 'n' if you don't want to create a new project, 'm' if you want to change/create only the Makefile, [another project name] if you want to create a new folder structure and makefile: " % project_path)
        if(answer=='n'):
            sys.exit('Project was not created')
        elif(answer=='m'):
            sys.stderr.write('Makefile was changed/added to the existing project %s\n' % project_path)
            break
        else:
            project_path = os.path.abspath(answer)
                    
                    

#######################################################################################################################
# Log project info
with open(os.path.join(project_path, 'log', 'info.txt'), 'w') as f:
    timestr = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    f.write("Time created: %s\n\nProject call: python %s\n" % (timestr, " ".join(sys.argv)) );



########################################################################################################################
## Main function to create one-sample Makefile
def makefile_local(m_input,  control):
    #print(m_input)
    #print(control)
    todel = []
    final_files = []
    mlist=[];
    if(type(m_input) == str):
        name = os.path.basename(m_input).split(".")[0];
    else:
        name = os.path.basename(m_input[0]).split(".")[0];
    
    log_dir = os.path.join('log', name)
    os.makedirs(os.path.join(project_path, log_dir), exist_ok=True);
    log_file = os.path.join(log_dir, "log.txt")
    

    # Processing bowite2 settings
    bs_list = get_bowtie_call(bowtie_settings, args.bowtie_args, args.index, m_input, name, threads=args.threads, reads_format=args.reads_format)
    bs_list = ['echo', '\'###bowtie\'', '>', log_file, ';'] + bs_list + ['2>> >(tee -a %s>&2)' % log_file]

    # Map reads with bowtie2
    input_files = m_input
    output_files = os.path.join('sam', '%s.sam' % name)
    todel.append(output_files)
    script = bs_list
    #print(script)
    mlist.append(dependence(input_files, output_files, script))
    
    # Convert mappings into coverage
    input_files = output_files;
    output_files = [os.path.join('coverage', '%s.%s.bed' % (name, x)) for x in ['plus', 'minus']]
    arguments = {'--genome': args.genome, '--outstat': log_dir, '--outcoverage': 'coverage', '--ambiguous': args.ambiguous}
    if(args.collapsed):
        arguments['--collapsed'] = 'True'
    script = get_script('get_sam_stat_paired.py', mapping_package, arguments=arguments, inp = input_files)
    mlist.append(dependence(input_files, output_files, script));   
    
    
    # UCSC coverage
    for of, strand in zip(output_files, ['plus', 'minus']):
        trackopts = "\'track name=%s_%s description=\"CHAP seq genomic coverage for sample %s\" %s\'" % (name, strand, name, " ".join(["%s=%s" % x for x in coverage_settings['trackopts'].items()]))
        input_files = of
        output_files = os.path.join('ucsc', '%s_%s.bedgraph' % (name, strand))
        final_files.append(output_files)
        script = get_script('coverage2bedgraph.py', mapping_package, arguments={'--multiplier': coverage_settings['multiplier'], '--convert': True, '--trackopts': trackopts}, inp = input_files, out = output_files)
        mlist.append(dependence(input_files, output_files, script));
    
    
    # Assign TPM and annotate the mappings
    
    input_files = [os.path.join('coverage', '%s.%s.bed' % (name, x)) for x in ['plus', 'minus']]
    output_files = os.path.join('transcripts', '%s.gff' % name)
    covpath = output_files
    script = get_script('assign_mappings.py', mapping_package, inp = input_files, out = output_files, arguments={'--transcripts': args.annotation, '--logdir': log_dir} )
    mlist.append(dependence(input_files, output_files, script));
            


    final_files.append(output_files)
    #Get header and cleaner for the makefile
    mlist.insert(0, get_header(final_files))
    todel = "\n\t".join( ['rm %s' % x  for x in todel] )
    mlist.append('clean:\n\t%s\n' % todel);

    return "\n\n".join(mlist), name, final_files


#######################################################################################################################
#Create Makefiles

if( os.path.isdir(args.reads[0]) ):
    input_list = []
    for path in args.reads:
        input_list.extend(list(sorted( [os.path.join(path, x) for x in listdir(path) if isfile(os.path.join(path, x))] )))
    input_list.sort()
else:
    input_list = args.reads;
if(args.paired):
    input_list = [(input_list[2*x], input_list[2*x+1]) for x in range(int(len(input_list)/2))]


    
    
if(args.control):
    if( os.path.isdir(args.control[0]) ):
        control_list = []
        for path in args.control:
            control_list.extend(list(sorted( [os.path.join(path, x) for x in listdir(path) if isfile(os.path.join(path, x))] )))
        control_list.sort()
    else:
        control_list = args.control;
    if(args.paired):
        control_list = [(control_list[2*x], control_list[2*x+1]) for x in range(int(len(control_list)/2))]
else:
    control_list = [None]*len(input_list);
        
            

if(args.paired):
    sample_names = [os.path.basename(x[0]).split(".")[0] for x in input_list]
else:
    sample_names = [os.path.basename(x).split(".")[0] for x in input_list]
mf_names = []
all_outputs = []
for m_input, control in zip(input_list, control_list):
    local_makefile, mname, local_output =  makefile_local(m_input, control)
    mname = 'makefile_%s' % mname
    mf_names.append(mname);
    all_outputs.append(local_output);
    with open(os.path.join(project_path, 'makefiles', mname), 'w') as mf:
        mf.write(local_makefile);



#######################################################################################################################  
### Create Master makefile ###


mlist = [];
final_files = []
mf_multipath = [os.path.join(project_path, 'makefiles', x) for x in mf_names] 

for mf_name, mf_path, input_names in zip(mf_names, mf_multipath, input_list):
    if(type(input_names) == str):
        input_names = [input_names]
    input_files = [mf_path] + list(input_names)
    output_files = mf_name;
    script = ["$(MAKE)", '-f', '$<']
    mlist.append(dependence(input_files, output_files, script))

if(args.multi):
    
    input_files = mf_names
    output_files = [os.path.join('log', 'report.%s' % x) for x in ('html', 'tsv')] 
    final_files.extend(output_files)
    script = get_script('log_html_total.py', rnaseq_package, arguments={'--css': os.path.join(html_lib, 'table.css'), '--js': os.path.join(html_lib, 'table.js'), '--name': args.name, '--order': sample_names, '--outdir': 'log'}, inp = 'log')
    mlist.append(dependence(input_files, output_files, script))    
    
   

        
        
#makefile header    
mlist.insert(0, get_header(final_files, phonyfiles=mf_names))
# makefie cleaner
mlist.append( 'clean:\n%s' %  ("\n".join(["\t$(MAKE) -f %s clean" % x for x in mf_multipath])) );
    

with open(os.path.join(project_path, 'Makefile'), 'w') as mf:
    mf.write("\n\n".join(mlist));


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


		
		
		
