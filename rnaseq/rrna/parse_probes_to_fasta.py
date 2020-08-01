import sys

with open(sys.argv[1]) as f:
    next(f);
    for n, l in enumerate(f, start=1):
        a = l.strip().split("\t")
        print(">probe_%d\n%s" % (n, a[1]))
