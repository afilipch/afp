'''Collection of functions supporting Makefile generation'''
import os
import sys

#chiflex_package = os.path.abspath(os.path.join(args.chiflex, 'chiflex'));
#splicing_package = os.path.abspath(os.path.join(args.chiflex, 'splicing'))
#interaction_package = os.path.abspath(os.path.join(args.chiflex, 'interaction'))
#mirna_package = os.path.abspath(os.path.join(args.chiflex, 'mirna'))

#chiflex_package = r'~/nrlbio/chiflex'


def get_script(script, package, arguments={}, inp = '', out = None, log=''):
    '''Example: get_script(something.py, {'--a': 7}, 'inp.txt', 'out.txt') will output: ('python [chiflex_package]/something.py'), 'inp.txt', '--a', '7', '>', 'out.txt')'''
    if(isinstance(inp, str)):
        input_files = inp
    else:
        input_files = " ".join(inp);

    l = [" ".join(('python', os.path.join(package, script) )), input_files]
    for k,v in arguments.items():
        if(v==False and v!=0):
            pass;
        elif(v==True and v!=1):
            l.append(k)
        elif(hasattr(v, '__iter__') and not isinstance(v, str)):
            l.append(k)
            l.append(" ".join([str(x) for x in v]))
        else:
            l.append(k)
            l.append(str(v))
    if(out):
        l.append(">")
        l.append(out)
    if(log):
        l.append('2>> >(tee -a %s>&2)' % log)
    return l
            


def dependence(input_files, output_files, script):
    '''Creates Makefile dependence'''
    if(isinstance(input_files, str)):
        inp = input_files
    else:
        inp = " ".join(input_files);
    if(isinstance(output_files, str)):
        out = output_files
    else:
        out = " ".join(output_files);
    return "%s: %s\n\t%s" % (out,inp, " ".join([str(x) for x in script]))



def get_header(output_files, phonyfiles=[]):
    if(isinstance(output_files, str)):
        ofs = output_files 
    else:
        ofs = " ".join(list(output_files))
    if(phonyfiles):    
        if(isinstance(phonyfiles, str)):
            phs = output_files 
        else:
            phs = " ".join(list(phonyfiles))

    if(phonyfiles):
        return "SHELL=/bin/bash\n.DELETE_ON_ERROR:\n\nall_clean: all clean\nall: %s\n.PHONY: %s all all_clean clean" % (" ".join((ofs, phs)), phs)
    else:
        return "SHELL=/bin/bash\n.DELETE_ON_ERROR:\n\nall: %s" % ofs



def get_bowtie_call(settings, arguments, reference, reads, project_name, threads=1, reads_format='U'):
    for bo in arguments:
        try:
            name, value = bo.split("=");
        except(ValueError):
            sys.exit("Bowtie options are provided in a malformatted way. Please see help for example\n") 
                
        if(name in settings):
            settings[name] = (value, settings[name][1])
        else:
            sys.stderr.write("provided option \'%s\' is currently not supported and will be ignored\n" % name) 
            
    settings['x'] = os.path.abspath(reference), '-'
    settings['S'] = os.path.join('sam', '%s.sam' % project_name), '-'
    if(isinstance(reads, str)):
        settings[reads_format] = reads, '-';
    else:
        settings['1'] = reads[0], '-'
        settings['2'] = reads[1], '-'
        

    settings['p'] = threads, '-';
    
    bs_list = ['bowtie2'];
    for k, v in settings.items():
        if(v[0]=='True'):
            bs_list.append(v[1] + k);
        elif(v[0]=='False'):
            pass;
        else:
            bs_list.append(v[1] + k);
            bs_list.append(v[0])
                    
    return bs_list;



def get_bowtie_help(bowtie_configurations):
    bowtie_help_list = [];
    for mode, settings in bowtie_configurations.items():
        bowtie_help_list.append("\n%s:" % mode)
        bowtie_help_list.append("[%s]" % " ".join(["%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in settings.items()]))
    return " ".join(bowtie_help_list)

