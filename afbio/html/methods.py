'''Collections of functions and classes for html programming in python3'''

def add_ucsc(interval, session, flank=0, chr_dict={}):
    if(chr_dict):
        return 'https://genome.ucsc.edu/s/%s?position=%s:%d-%d' % (session, chr_dict[interval.chrom], interval.start-flank, interval.stop+flank)
    else:
        return 'https://genome.ucsc.edu/s/%s?position=chr1:%d-%d' % (session, interval.start-flank, interval.stop+flank)

