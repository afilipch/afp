import os
import sys


with open(sys.argv[1]) as f:
    header = next(f)
    for l in f:
        a = l.strip().split("\t")
        print(">%s:%s-%s\n%s" % tuple(a[:4]))
