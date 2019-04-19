#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Creates makefile and directory structure for chip/chap seq analysis'''


import argparse
import sys
import os

from afbio.config.config import load_config;
from afbio.makefiles import dependence, get_header, get_bowtie_call, get_script, get_bowtie_help

#Load default parameters for bowtie2
conf = load_config('chipchap')

#Bowtie options preliminary parsing
bowtie_settings = conf['bowtie'];
peak_detection_settings = conf['peak_detection']
coverage_settings = conf['coverage']
annotation_settings = conf['annotation']
region_settings = conf['region']

bowtie_help_str = "[%s]" % " ".join(["%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in bowtie_settings.items()])


parser = argparse.ArgumentParser(description='Creates makefile and directories for chiflex project')#, formatter_class = argparse.RawTextHelpFormatter);
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = os.path.abspath, help = "Path to the project folder. If folder does not exist, it will be created");
parser.add_argument('--package', nargs = '?', type = os.path.abspath, required = True, help = "Path to the afp package");
parser.add_argument('--reads', nargs = '+', type = os.path.abspath, help = "Path to sequencing reads. fastq/fasta file. Paired-end reads must be provided consecutively");
parser.add_argument('--index', nargs = '?', type = os.path.abspath, help = "Path to the mapping reference bowtie2 indices");
parser.add_argument('--genome', nargs = '?', type = os.path.abspath, help = "Path to the mapping references in fasta format");

#parser.add_argument('--reads', nargs = '+', type = os.path.abspath, required = True, help = "Path to sequencing reads. fastq/fasta file. Paired-end reads must be provided consecutively");
#parser.add_argument('--index', nargs = '?', type = os.path.abspath, required = True, help = "Path to the mapping reference bowtie2 indices");
#parser.add_argument('--genome', nargs = '?', type = os.path.abspath, required = True, help = "Path to the mapping references in fasta format. If set, the sequences will be assigned to resulting interactions as well as energy and binding pattern(pattern of paired nucleotides on a left chimeric part)");

parser.add_argument('--coverage', nargs = '+', type = os.path.abspath, help = "Path to the already generated coverage files. If provided, --reads argument is obsolete and ignored");

#Options for the mapping result postprocessing
parser.add_argument('--collapsed', nargs = '?', default = False, const=True, type = bool, help = "If set, reads are assumed collapsed with collapse.pl script. Read count appendix of the read id will be used to calculate read support of the interactions")
parser.add_argument('--paired', nargs = '?', default = False, const=True, type = bool, help = "If set, reads are assumed to be paired-end")

#Options for the output control
parser.add_argument('--only_makefile', nargs = '?', default = False, const = True, type = bool, help = "If set, a new makefile is created, but not folders");

#Options controlling postprocessing
parser.add_argument('--multi', nargs = '?', default = False, const = True, type = bool, help = "If set, consensus regions will be generated and explored");
parser.add_argument('--annotation', nargs = '?', type = os.path.abspath, help = "Path to genbank annotation file in gff format.");

#bowtie2 options
parser.add_argument('--bowtie_args', nargs = '+', default = [], type = str, help = "Bowtie settings for the first round of mapping. For example, if one wants to set \'-p 4\', use \'--local\' alignment mode, but not \'--norc\' option then \'p=4 local=True norc=False\' should be provided. Given attributes replace default(for Chiflex, NOT for Bowtie) ones. Default settings for the modes are: %s" % bowtie_help_str)
args = parser.parse_args();


#######################################################################################################################
# Process input options


if(args.reads):
    if(not (args.genome and args.index) ):
        parser.error('If --reads are provided, --genome and --index must be set as well');
elif(not args.coverage):
    parser.error('Either --reads or --coverage argument must be provided');
    



#######################################################################################################################
# Set paths' constants

chap_package = os.path.join(args.package, 'chap')
mapping_package = os.path.join(args.package, 'mapping')


#######################################################################################################################
# Create folders

project_path = args.path
folders = ['sam', 'coverage', 'statistics', 'log', 'peaks', 'makefiles', 'regions', 'ucsc']


while(not args.only_makefile):
    try:
        os.makedirs(project_path);
        for folder in folders:
            os.mkdir(os.path.join(project_path, folder));
        break;
    except:
        answer = input("\nProject directory \'%s\' is currently exists, please type 'N' if you don't want to create a new project, 'MO' if you want to change/create only the Makefile, [another project name] if you want to create a new folder structure and makefile: " % project_path)
        if(answer=='N'):
            sys.exit('Project was not created')
        elif(answer=='MO'):
            sys.stderr.write('Makefile was changed/added to the existing project %s\n' % project_path)
            break
        else:
            project_path = os.path.abspath(answer)
                    
                    




########################################################################################################################
## Main function to create one-sample Makefile
def makefile_local(m_input, coverage_mode, multi=False):
    mlist=[];
    if(type(m_input) == str):
        name = os.path.basename(m_input).split(".")[0];
    else:
        name = os.path.basename(m_input[0]).split(".")[0];

    if(not coverage_mode):
        # Processing of the left chimeric part bowite2 settings
        bs_list = get_bowtie_call(bowtie_settings, args.bowtie_args, args.index, m_input, name)

        # Map reads with bowtie2
        input_files = m_input
        output_files = os.path.join('sam', '%s.sam' % name)
        script = bs_list
        mlist.append(dependence(input_files, output_files, script))
        
        # Convert mappings into coverage
        input_files = output_files;
        output_files = [os.path.join('coverage', '%s.%s.bed' % (name, x)) for x in ['minus', 'plus']]
        script = get_script('get_sam_stat_paired.py', mapping_package, arguments={'--genome': args.genome, '--outstat': 'statistics', '--outcoverage': 'coverage'}, inp = input_files)
        mlist.append(dependence(input_files, output_files, script));   
        
        # Merge coverages coming from different strands
        input_files = output_files;
        output_files = os.path.join('coverage', '%s.bed' % name)
        covpath = output_files
        script = get_script('merge_coverages.py', chap_package, inp = input_files, out = output_files)
        mlist.append(dependence(input_files, output_files, script));
    else:
        covpath = m_input;
        output_files = m_input
    
    # Detect peaks
    input_files = output_files
    output_files = [os.path.join('peaks', '%s.raw.bed' % name), os.path.join('log', '%s.convolution.bed' % name), os.path.join('log', '%s.convolution.png' % name)]
    script = get_script('detect_peaks.py', chap_package, arguments={'--widthfactor': peak_detection_settings['widthfactor'], '--meanmult': peak_detection_settings['meanmult'], '--convolution': output_files[1] , '--plot': output_files[2]}, inp = input_files, out = output_files[0])
    mlist.append(dependence(input_files, output_files, script));  
    
    # Filter peaks
    input_files = output_files[0]
    output_files = [os.path.join('peaks', '%s.filtered.bed' % name), os.path.join('log', '%s.assigned.tsv' % name), os.path.join('statistics', '%s.filtering.png' % name)]
    filtered_path = output_files[0]
    script = get_script('filter_peaks.py', chap_package, arguments={'--zscore': peak_detection_settings['zscore'], '--minmedian': peak_detection_settings['minmedian'], '-ap': output_files[1] , '--plot': output_files[2], '--coverage': covpath}, inp = input_files, out = output_files[0])
    mlist.append(dependence(input_files, output_files, script));
    
    # Normalize coverage
    input_files = [output_files[1], covpath]
    output_files = os.path.join('coverage', '%s.normalized.bed' % name)
    normed_covpath = output_files;
    script = get_script('normalize_coverage.py', chap_package, arguments={'--zscore': coverage_settings['zscore'], '--mode': coverage_settings['mode'], '--coverage': input_files[1]}, inp = input_files[0], out = output_files)
    mlist.append(dependence(input_files, output_files, script));
    
    # Create UCSC tracks
    # python ~/afp/mapping/coverage2bedgraph.py coverage/time10_neb.normalized.bed --trackopts 'track name=time10_neb description="CHAP coverage for time10_neb sample" color=60,150,237 type=bedGraph visibility=2' --multiplier 100 --convert > todel.bedgrap
    trackopts = "\'track name=%s description=\"CHAP seq genomic coverage for sample %s\" %s\'" % (name, name, " ".join(["%s=%s" % x for x in coverage_settings['trackopts'].items()]))
    input_files = output_files
    output_files = os.path.join('ucsc', '%s.bedgraph' % name)
    script = get_script('coverage2bedgraph.py', mapping_package, arguments={'--multiplier': coverage_settings['multiplier'], '--convert': True, '--trackopts': trackopts}, inp = input_files, out = output_files)
    mlist.append(dependence(input_files, output_files, script));
    
    
    
    # Annotate peaks
    if(not multi):
        input_files = [filtered_path, output_files]
        output_files = os.path.join('peaks', '%s.annotated.gff' % name)
        script = get_script('annotate.py', chap_package, arguments={'--maxshift': annotation_settings['maxshift'], '--flen': annotation_settings['flen'], '--coverage': input_files[1], '--genes': args.annotation}, inp = input_files[0], out = output_files)
        mlist.append(dependence(input_files, output_files, script)); 
    else:
        output_files = [output_files] + [filtered_path]
    
    #Get header and cleaner for the makefile
    mlist.insert(0, get_header(output_files))
    mlist.append('clean:\n\techo "nothing to clean."\n');

    return "\n\n".join(mlist), name, [normed_covpath] + [filtered_path]


#######################################################################################################################
#Create Makefiles

if(args.reads):
    if(args.paired):
        input_list = [(args.reads[2*x], args.reads[2*x+1]) for x in range(int(len(args.reads)/2))]
    else:
        input_list = args.reads;
    coverage_mode = False
else:
    input_list = args.coverage;
    coverage_mode = True

#input_list = [x for x in input_list]
#sys.stderr.write("%s\n" % input_list)
sample_names = [os.path.basename(x[0]).split(".") for x in input_list]



 
mf_names = []
all_outputs = []
for m_input in input_list:
    local_makefile, mname, local_output =  makefile_local(m_input, coverage_mode, args.multi)
    mname = 'makefile_%s' % mname
    mf_names.append(mname);
    all_outputs.append(local_output);
    with open(os.path.join(project_path, 'makefiles', mname), 'w') as mf:
        mf.write(local_makefile);

#######################################################################################################################  
# Create Master makefile


mlist = [];
global_output = []
mf_multipath = [os.path.join(project_path, 'makefiles', x) for x in mf_names] 

for mf_name, mf_path, input_names in zip(mf_names, mf_multipath, input_list):
    if(type(input_names) == str):
        input_names = [input_names]
    input_files = [mf_path] + list(input_names)
    output_files = mf_name;
    script = ["$(MAKE)", '-f', '$<']
    mlist.append(dependence(input_files, output_files, script))

if(args.multi):
    all_coverages = [x[0] for x in all_outputs]
    all_peaks = [x[1] for x in all_outputs]
    
    input_files = mf_names
    output_files = os.path.join('regions', 'regions.gff')
    raw_regions = output_files
    script = get_script('merge_peaks.py', chap_package, arguments={'--coverage': all_coverages, '--zscore': region_settings['zscore'], '--flank': region_settings['flank']}, inp = all_peaks, out = output_files)
    mlist.append(dependence(input_files, output_files, script));
    
    input_files = output_files
    output_files = os.path.join('statistics', 'peaks.correlation.png')
    global_output.append(output_files)
    script = get_script('correlate_peaks.py', chap_package, arguments={'--min-zscore': region_settings['min-zscore'], '--names': sample_names, '--plot': output_files}, inp = input_files)
    mlist.append(dependence(input_files, output_files, script));
    
    #python /home/a_filipchyk/afp/chap/annotate.py regions/regions.gff --maxshift 50 --flen 50  --genes /home/a_filipchyk/genomic_data/coryne/annotation/improved_annotation_2017.gff
    if(args.annotation):
        input_files = raw_regions
        output_files = os.path.join('regions', 'regions.annotated.gff')
        script = get_script('annotate.py', chap_package, arguments={'--maxshift': region_settings['maxshift'], '--flen': region_settings['flank'], '--genes': args.annotation}, inp = input_files, out = output_files)
        mlist.append(dependence(input_files, output_files, script));    
    
    if(type(output_files) == str):
        global_output.append(output_files);
    else:
        global_output.extend(output_files);
        
        
#makefile header    
mlist.insert(0, get_header(global_output, phonyfiles=mf_names))
# makefie cleaner
mlist.append( 'clean:\n%s' %  ("\n".join(["\t$(MAKE) -f %s clean" % x for x in mf_multipath])) );
    

with open(os.path.join(project_path, 'Makefile'), 'w') as mf:
    mf.write("\n\n".join(mlist));


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


		
		
		
