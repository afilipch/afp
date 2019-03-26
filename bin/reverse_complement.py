#! /home/a_filipchyk/soft/home/a_filipchyk/anaconda3/bin/python
'''Outputs reverse complement for the provided sequences'''

import pybedtools
import sys
from Bio.Seq import reverse_complement


for seq in sys.argv[1:]:
    print(reverse_complement(seq));