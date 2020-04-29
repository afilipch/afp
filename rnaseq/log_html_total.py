import sys
import os
import dominate
import argparse
from collections import defaultdict
from dominate.tags import *
from dominate.util import raw


from afbio.html.methods import add_ucsc
from math import log

parser = argparse.ArgumentParser(description='Generates html report of CHAP analyses by chipchap.py');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the root log folder");
parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
parser.add_argument('--name', nargs = '?', default="experiment", type = str, help = "Name of the experiment");
parser.add_argument('--order', nargs = '+', required=True, type = str, help = "Names of the provided samples in the order to be reported");
#parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
args = parser.parse_args();

### SET CONSTANTS ###

bowtie_labels = ["Sample", "Total Reads", "Mapped Uniq[%]", "Mapped Ambiguously[%]"] #, "Mapped non-uniquely", "Unmapped reads", "Mapped Discordantly"]
bowtie_dtypes = [0, 1, 1, 1]#, 1, 1, 1]
types_labels = ['protein coding', 'rRNA', 'intergenic', 'tRNA', 'pseudogene', 'misc RNA', 'rRNA species']

labels = bowtie_labels + ["%s [%%]" % x for x in types_labels]
dtypes = bowtie_dtypes + [1]*len(types_labels)


### FUNCTIONS ###
def process_log(log_file):
    section2text = {} #defaultdict(list);
    section = ""
    with open(log_file) as f:
        for l in f:
            if(l.startswith("###")):
                section = l.strip()[3:];
                section2text[section] = [];
            elif(section):
                section2text[section].append(l.strip())
    return section2text;


def get_bowtie(section2text, btype):
    bowtie_list = [];
    ltemp = section2text[btype][1:8]
    for l in (ltemp[:4] + ltemp[6:]):
        bowtie_list.append(int(l.split(" ")[0]))
    if( len(bowtie_list)>4 ):   
        bowtie_list[1] -= bowtie_list[4]
    bowtie_list[1], bowtie_list[2], bowtie_list[3] = bowtie_list[2], bowtie_list[3], bowtie_list[1]
    return bowtie_list


def get_types(section2text):
    filter_dict = {};
    rrna = []
    for l in (section2text['types']):
        if(l.startswith('#')):
            rrna.append("%s=%s" % tuple(l[1:].split("\t")) );
        else:
            key, value = l.split("\t")
            key = key.replace("_", " ")
            filter_dict[key] = value;
    filter_dict['rRNA species'] = "; ".join(rrna)
    return filter_dict
    

def process_folder(folder):
    data_dict = {}
    section2text = process_log(os.path.join(folder, "log.txt"))
    #print(section2text)
    data_dict['bowtie'] = get_bowtie(section2text, 'bowtie')
    if('bowtie_control' in section2text):
        data_dict['bowtie_control'] = get_bowtie(section2text, 'bowtie_control')
    data_dict['types'] = get_types(section2text)
    return data_dict



#print([x for x in os.walk(args.path)])

folders = [os.path.join(args.path, x) for x in args.order]
data = [process_folder(x) for x in folders]
names = [os.path.basename(x) for x in folders]






#sys.exit()





        






########### HTML SECTION ###########
doc = dominate.document(title="%s: CHAP analyses overview" % args.name)
with open(args.css) as f:
    _style = f.read()
with doc.head:
    style(_style)
    
with open(args.js) as f:
    plain_script = f.read()
    
with doc:
    p(strong("Mapping Results"))
    #input(type="text", id="myInput", onkeyup="my_search(4)", placeholder="Filter relation to CgpS..")
    #input(type="number", id="myInputGreater", onkeyup="my_filter_greater(6)", placeholder="Filter AT content greater than..")
    with table(id = "myTable") as _table:
        _tr = tr()
        _tr.add([ th(x[1][0], onclick='sortTable(%d, %d)' % (x[0], x[1][1])) for x in enumerate(zip(labels, dtypes))  ])
        for name, sample_data in zip(names, data):
            #print(sample_data['bowtie'])
            with tr():
                td(raw('<a href=%s target="_blank">%s</a>' % (os.path.join(name, "report.html"), name) )) 
                td(sample_data['bowtie'][0])
                td("%1.1f" % (sample_data['bowtie'][1]/sample_data['bowtie'][0]*100))
                td("%1.1f" % (sample_data['bowtie'][2]/sample_data['bowtie'][0]*100))
                for label in types_labels:
                    td(sample_data['types'][label])




                
    _script = script(raw(plain_script), type='text/javascript')

print(doc.render());





