'''Collections of functions and classes for html programming in python3'''

def add_ucsc(interval, session, flank=0):
    ucsc_link = 'https://genome.ucsc.edu/s/%s?position=chr1:%d-%d' % (session, interval.start-flank, interval.stop+flank)
    return ucsc_link
