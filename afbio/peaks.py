import numpy as np;
from afbio.sequencetools import sliding_window
from multiprocessing.dummy import Pool
import sys;


        
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


def convolute(arr, kernel, peakwidth, threads = 4, takepeak=True):
    scaled = kernel2scale(kernel, peakwidth, takepeak);
    def _local_convolution(window):
        return sum(window*scaled) 

    
    wsize = len(scaled)
    tail = np.zeros((wsize-1)//2)
    extended = np.concatenate((tail, arr, tail))
    pool = Pool(threads)
    #windows = [np.array(x) for x in sliding_window(extended, wsize)]
    #sys.stderr.write('f\n')
    
    return [_local_convolution(np.array(x)) for x in sliding_window(extended, wsize)]
    #return pool.imap(_local_convolution, sliding_window(extended, wsize))


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
    
    
            
            
