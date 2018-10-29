import math;
import os;
import sys;
import numpy as np;
from collections import Counter;
from itertools import islice

def sliding_window(seq, n):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def gauss_second_kernel(sigma=1.0, truncate=4.0, steps=100000):
    res = [];
    c = 1.0/(math.sqrt(2*math.pi)*sigma**5)
    for x in np.linspace(-sigma*truncate, sigma*truncate, num = steps):
        res.append(math.exp(-x**2/(2*sigma**2))*(sigma**2-x**2)*c);
    return res;

def kernel2scale(kernel, peakwidth):
    bw = 0;
    km12 = max(kernel)/2.0;
    for c, el in enumerate(kernel):
        if(el>km12):
            bw = (0.5 - float(c)/len(kernel))*2;
            break;
    flank_size = int(round(peakwidth/(bw*2)));
    step = len(kernel)/(flank_size*2+1)
    ends = [int(round(x)) for x in np.cumsum([step]*(flank_size*2+1))]
    starts = [0] + ends[:-1]
    
    scaled = [];
    for s,e in zip(starts, ends):
        a = kernel[s:e];
        val = np.mean(a)
        #pos = int((e-s)/2)
        #tans = [0]*(e-s)
        #tans[pos] = val
        scaled.append(val)
    #ans = [];
    #curarr = [];
    #for c, val in enumerate(kernel
    
    return np.array(scaled);
            
    

def convolute(arr, kernel, peakwidth):
    scaled = kernel2scale(kernel, peakwidth);
    wsize = len(scaled)
    tail = np.zeros((wsize-1)//2)
    extended = np.concatenate((tail, arr, tail))
        #signal = [];
    for c, window in enumerate(sliding_window(extended, wsize)):
        yield sum(np.array(window)*scaled)
    #return signal


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


###Aliases
ndg = gauss_second_kernel

        
   
   
   
#TESTING    
if(__name__ == "__main__"):       

    mu_arr = [200, 260, 550, 800, 970]; 
    sigma_arr = [20, 25, 15, 30, 20];
    arrs = [];
    for mu, sigma in zip(mu_arr, sigma_arr):
        arrs.append(np.random.normal(loc=mu, scale=sigma, size = 100000))
        
    arr = np.array([int(round(x)) for x in np.concatenate(arrs)])
    arrcounter = Counter(arr)
    arrcounter = [arrcounter.get(x,0) for x in range(max(arrcounter.keys()) + 1)]
    arrcounter.extend([0]*100)
    arrcounter = np.array(arrcounter)
    #print(arrcounter)
    #print(len(arr_filtered))


    kernel = gauss_second_kernel(1);        
    #print(kernel)    

    #scaled_kernel = kernel2scale(kernel, 20)
    signal = convolute(arrcounter, kernel, 25)
    peaks = detect_peaks(signal)
    px, py = zip(*sorted(peaks.items()))
    #sys.exit()

    #for w in sliding_window('123456789',3):
        #print(w)

    #PLOT    
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots()

    ax1.plot(signal, 'b-')
    ax1.set_xlabel("postion (nt)")
    ax1.set_ylabel('coverage', color='b')
    ax1.tick_params('y', colors='b')

    #ax2 = ax1.twinx()
    ax1.plot(px, py, 'r*')
    #ax2.set_ylabel('detected peaks', color='r')
    #ax2.tick_params('y', colors='r')

    fig.tight_layout()
    plt.show()
