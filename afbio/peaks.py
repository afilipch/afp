import sys
from collections import defaultdict
from multiprocessing import Pool

from pybedtools import BedTool
import numpy as np;

from afbio.sequencetools import sliding_window




        
def kernel2scale(kernel, peakwidth, takepeak):
    if(takepeak):
        bw = 0;
        km12 = max(kernel)/2.0;
        #sys.stderr.write("bu" + "\n")
        for c, el in enumerate(kernel):
            if(el>km12):
                bw = (0.5 - float(c)/len(kernel))*2;
                #sys.stderr.write("%1.1f\n" % bw )
                break;
        #sys.stderr.write("%1.1f\n" % peakwidth )
        flank_size = int(round(peakwidth/(bw*2)));
    else:
        flank_size = peakwidth//2;
    step = len(kernel)/(flank_size*2+1)
    ends = [int(round(x)) for x in np.cumsum([step]*(flank_size*2+1))]
    starts = [0] + ends[:-1]
    
    scaled = [];
    for s,e in zip(starts, ends):
        a = kernel[s:e];
        val = np.mean(a)
        scaled.append(val)
    
    return np.array(scaled);


def local_convolution(alist):
    return sum(alist[0]*alist[1]) 

def local_generator(extended, wsize, scaled):
    for window in sliding_window(extended, wsize):
        yield (window, scaled);

def convolute(arr, kernel, peakwidth, threads, takepeak=True):
    scaled = kernel2scale(kernel, peakwidth, takepeak);
    
    wsize = len(scaled)
    tail = np.zeros((wsize-1)//2)
    extended = np.concatenate((tail, arr, tail))
    
    result = [];
    with Pool(threads) as p:
        for pos, c in enumerate(p.imap(local_convolution, local_generator(extended, wsize, scaled), chunksize = threads*10)):
            #if(pos  and pos % 100000 == 0):
                #sys.stderr.write("%d nucleotides are already processed\n" % pos)
            result.append(c)
    return result


def detect_peaks(signal):
    locs = [];
    ispeak = False;
    start = 0;
    
    for c, window in enumerate(sliding_window(signal, 3)):
        if((window[0] > window[1] < window[2])):
            if(ispeak):
                end = c+1;
                ispeak = False;
                locs.append((start, top, end, height))
                start = c+1;
            else:
                start = c+1;
                
        if(window[1]*window[2] < 0):
            if(ispeak):
                end = c+2;
                ispeak = False;
                locs.append((start, top, end, height))
            else:
                start = c+2
                
        if((window[0] < window[1] > window[2]) and window[1]>0):
            top = c+1;
            ispeak = True;
            height =  window[1]
    return locs;


def findpeak(region, threshold):
    mh = max(region)
    if(mh<threshold):
        return None;
    
    else:
        mp = np.argmax(region)
        half = mh/2.0
        
        for p, el in enumerate(region[mp+1:]):
            if(el>half):
                end = mp + 1 + p
                break;
        else:
            return None

        for p, el in enumerate(region[mp-1::-1]):
            if(el<half):
                start = mp - 1 - p
                break;
        else:
            return None
        
        return start, mp, end;
    
def estimate_bandwidth(arr_coverage, meanmult):
    peak_threshold = np.mean(arr_coverage)*meanmult;
    toppeaks = set();
    steps = [880, 1430, 2000];
    for step in steps:
        for k in range(len(arr_coverage)//step):
            region = arr_coverage[k*step:(k+1)*step]
            peak = findpeak(region, peak_threshold)
            if(peak):
                apeak = tuple([x+k*step for x in peak]);
                toppeaks.add(apeak)
    lengths = [x[2] - x[0] for x in toppeaks]
    return np.mean(lengths), len(toppeaks)


def recenter_based_on_coverage(peak, coverage, fold):
    mpos = np.argmax(coverage);
    threshold = max(coverage)*fold;
    for c, val in enumerate(coverage[mpos:]):
        if(val<=threshold):
            right = c+mpos;
            break;
    else:
        right = len(coverage)-1
        
    for c, val in enumerate(coverage[mpos::-1]):
        if(val<=threshold):
            left = mpos-c;
            break;
    else:
        left = 0
            
    return (left, mpos, right);





def _find_shared_peaks_chromosome(bedtools_chr, maxd):
    res = [];
    stat_counts = []
    marked = [];
    for c, intervals in enumerate(bedtools_chr):
        for interval in intervals:
            marked.append((c, interval))
    marked.sort(key = lambda x: int(x[1].name));
    selections = [];
    current_selection = [marked[0]];
    for m in marked[1:]:
        if(current_selection and int(m[1].name) - int(current_selection[0][1].name) <= maxd):
            current_nums = [x[0] for x in current_selection]
            if(m[0] not in current_nums):
                current_selection.append(m)
        else:
            res.append(current_selection)
            stat_counts.append(len(current_selection))
            current_selection = [m]
    else:
        res.append(current_selection)
        stat_counts.append(len(current_selection))
        
    return res, stat_counts



def find_shared_peaks(multipath, maxd):
    #Split by chromosome
    
    chr2bedtools = defaultdict(list);
    for intervals in [BedTool(x) for x in multipath]:

        temp_d = defaultdict(list);
        for interval in intervals:
            temp_d[interval.chrom].append((interval));
        for chrom, local_intervals in temp_d.items():
            chr2bedtools[chrom].append(local_intervals)
            
    bedtools_list = [x[1] for x in sorted(chr2bedtools.items(), key = lambda x: x[0])]

    stat_total_counts = []
    res_total = [];
    for bedtools_chr in bedtools_list:
        res, stat_counts = _find_shared_peaks_chromosome(bedtools_chr, maxd)
        res_total.extend(res)
        stat_total_counts.extend(stat_counts);
        
    return res_total, stat_total_counts


def shared_peaks_stat_to_string(stat_total_counts, size):
    l = []
    l.append("number of peaks per merged\tnumber of merged peaks\tfraction [%]\n")    
    for s in range(1, size+1):
        count = stat_total_counts.count(s);
        l.append("%d\t%d\t%1.1f\n" % (s, count, count/len(stat_total_counts)*100))
    return "".join(l)
    
    
    
    
    
    
    
    
        
        
    
    
    
            
            
