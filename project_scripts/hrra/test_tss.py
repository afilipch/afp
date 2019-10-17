import sys
from pybedtools import BedTool

b1 = BedTool(sys.argv[1])
b2 = dict([ (x.name, x) for x in BedTool(sys.argv[2]) ])

for p1 in b1:
    p2 = b2.get(p1.name)
    if(not p2):
        print(p1);
        continue
    if(p1.attrs['gene'] != p2.attrs['gene']):
        print(p1, p2)

