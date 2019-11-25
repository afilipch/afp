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
#parser.add_argument('--genome', nargs = '?', required=True, type = str, help = "Path to the genome, fasta format");
parser.add_argument('--js', nargs = '?', required=True, type = str, help = "Path to javascript functions");
parser.add_argument('--css', nargs = '?', required=True, type = str, help = "Path to css style sheet");
#parser.add_argument('--ucsc', nargs = '?', required=True, type = str, help = "Name of the UCSC session");
#parser.add_argument('--flank', nargs = '?', default=60, type = int, help = "Peak plank length");
args = parser.parse_args();

### SET CONSTANTS ###

bowtie_labels = ["Sample", "Total Reads", "Mapped [%]"] #, "Mapped non-uniquely", "Unmapped reads", "Mapped Discordantly"]
bowtie_dtypes = [0, 1, 1]#, 1, 1, 1]
filter_labels = ["Total Peaks", "Score Threshold", "Filtered/Score", "Filtered/Coverage", "Re-centered peaks [%]"]
filter_dtypes = [1, 1, 1, 1, 1]

labels = bowtie_labels + filter_labels
dtypes = bowtie_dtypes + filter_dtypes


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
        
    bowtie_list[1] -= bowtie_list[4]
    bowtie_list[1], bowtie_list[2], bowtie_list[3] = bowtie_list[2], bowtie_list[3], bowtie_list[1]
    return bowtie_list

### Filter peaks processing
#coverage_list = [];

def get_filter(section2text):
    filter_list = [];
    for l in (section2text['filter_peaks']):
        if(l):
            filter_list.append(float(l.split("\t")[1]))
    return [ int(filter_list[6]), filter_list[5], int(filter_list[7]), int(filter_list[12]), filter_list[13] ]
    
    


def process_folder(folder):
    data_dict = {}
    section2text = process_log(os.path.join(folder, "log.txt"))
    data_dict['bowtie'] = get_bowtie(section2text, 'bowtie')
    data_dict['bowtie_control'] = get_bowtie(section2text, 'bowtie_control')
    data_dict["filter"] = get_filter(section2text)
    return data_dict



#print([x for x in os.walk(args.path)])

folders = [os.path.join(args.path, x) for x in os.listdir(args.path) if os.path.isdir(os.path.join(args.path, x))]
data = [process_folder(x) for x in folders]
names = [os.path.basename(x) for x in folders]






#sys.exit()





        






########### HTML SECTION ###########
doc = dominate.document(title="CHAP analyses overview for all the samples")
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
                td(name)
                td(sample_data['bowtie'][0])
                td("%1.1f" % (sample_data['bowtie'][1]/sample_data['bowtie'][0]*100))
                td(sample_data['filter'][0])
                td("%1.1f" % sample_data['filter'][1])
                td(sample_data['filter'][2])
                td(sample_data['filter'][3])
                td("%1.1f" % (sample_data['filter'][4]/sample_data['filter'][3]*100))
                #td(sample_data['bowtie'][0][2])
                #td(sample_data['bowtie'][0][3])
                #td(sample_data['bowtie'][0][4])

                
                
    _script = script(raw(plain_script), type='text/javascript')
    #_script.add_raw_string(plain_script)


print(doc.render());










'''357955 reads; of these:
  357955 (100.00%) were paired; of these:
    31369 (8.76%) aligned concordantly 0 times
    321179 (89.73%) aligned concordantly exactly 1 time
    5407 (1.51%) aligned concordantly >1 times
    ----
    31369 pairs aligned concordantly 0 times; of these:
      1049 (3.34%) aligned discordantly 1 time
91.53% overall alignment rate
'''
