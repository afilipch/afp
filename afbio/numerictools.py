# /usr/bin/python
'''collections of classes and functions to solve diverse numerical problems'''
from math import log
from itertools import chain, islice
import random

import numpy as np;
from bisect import bisect_left



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


def select_by_probability(items, probabilities):
    '''Randomly select an item accroding to the provided probabilities
            items iterable: items to select from
            probabilities iterable: probabilities of the items to be selected
    Returns object: randomly selected item
    '''   
    norm = sum(probabilities)
    probs = np.cumsum([x/norm for x in probabilities])
    r = random.random()
    return items[bisect_left(probs, r)]
    
    


#def cdf(iterable):
    #from collections import Counter
    #'''Constructs cdf from given iterable with numeric entries
            #iterable iterable: entries have to be of numeric type
    #Returns list: each element is a tuple of two elements. 1st is any value from iterable, 2nd number of entries in iterable, which are less or equal to the 1st element.
    #'''
    #l = [];
    #sum = 0;
    #for k, v in sorted(Counter(iterable).items(), key = lambda x: x[0]):
        #sum+=v;
        #l.append((k, sum));
    #return l;



def maxes(list_, key_function=lambda x: x):
    '''Returns all maximum(according to key function) elements in a given iterable.
    
        list_ sequence type: any sequence type to look for maximum elements in
        key_function function returns float: function to be called on the element of an iterable. The bigger the value function returns, the 'bigger' the elements
        
        Return tuple: 1st element is a list of maximum elements. 2nd element is an integer output of key_function value associated with maximum element.
    '''
    if(not list_):
        return [], None;
            
    m, max_list = key_function(list_[0]), []
    for s in list_:
        k = key_function(s)
        if k > m:
            m, max_list = k, [s]
        elif k == m:
            max_list.append(s)
    return max_list, m


def overlap(i1, i2):
    '''Calculates overlap of two intervals, maybe negative

        i1 iterable: 1d interval. 2-element iterable. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive)
        i2 iterable: another 1d interval. 2-element iterable. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive)
        
        Returns tuple: Overlap of two intervals, tuple:
            First element is start of interval(0-based inclusive). 
            Second element is end of interval(0-based exclusive).
    '''
    start = max(i1[0], i2[0])
    end = min(i1[1], i2[1])
    return start, end;


def distance(i1, i2):
    '''Calculates the distance between two intervals
    
        i1 iterable: 1d interval. 2-element iterable. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive)
        i2 iterable: another 1d interval. 2-element iterable. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive)
        
        Returns int:  distance between two intervals
    '''
    return max(i1[0], i2[0]) - min(i1[1], i2[1])
	
	
def merge_intervals(intervals, distance=0, assume_sorted=False):
    '''Yields lists of overlapping intervals
    
        intervals iterable: element is 2-element tuple(First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive))
        distance int: minimum overlap(max gap if negative) required
        assume_sorted bool: if True, "intervals" argument is treated as sorted(according to start position) iterable
        
        Yields list: list of overlapping intervals (2-element tuple. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive))
    '''
    if(hasattr(intervals, 'next')):
        first = intervals.next();
        l = intervals
            
    else:	
        l = list(intervals)
        if(not assume_sorted):
            l.sort(key = lambda x: x[0]);	
        first = l.pop(0);
    
    merged = [first];
    start, end = first;
    for i in l:
        s, e = overlap(i, (start, end))
        if(e - s >= distance):
            merged.append(i);
            end = max(i[1], end)
        else:
            yield merged;
            start, end = i;
            merged = [i];
    yield merged
		
		
def interval_intersection(i1, i2):
    '''Returns overlap of two intervals, only positive, or None
            
        i1 list|tuple: 1st element is the start of the interval(0-based), 2nd element is the end of the interval(exclusive)
        i2 list|tuple: 1st element is the start of the interval(0-based), 2nd element is the end of the interval(exclusive)
        
        Returns tuple|None: overlap of two intervals. 1st element start of the overlap(0-based), 2nd element end of the overlap(exclusive). None if there is no overlap
    '''	
    start, end = overlap(i1, i2);
    if(end>start):
        return start, end;
    else:
        return None;
            
		
def overlap_hyperrectangles(hr1, hr2):
    '''Returns intersection of hyperrectangles
        
        hr1 list: hyperrectangle. List of intervals(2-element iterables: first element is start of interval(0-based inclusive), second element is end of interval(0-based exclusive))
        hr2 list: another hyperrectangle. List of intervals(2-element iterables: first element is start of interval(0-based inclusive), second element is end of interval(0-based exclusive))
        
        Returns list: overlap of two hyperrectangles. List of intervals(2-element iterables: first element is start of interval(0-based inclusive), second element is end of interval(0-based exclusive))
    '''
    c = [];
    for i1, i2 in zip(hr1, hr2):
        o = interval_intersection(i1, i2);
        if(o):
            c.append(list(o));
        else:
            return None
    return c;
		
	
		
def dict2entropy(d):
    '''Calculates entropy for a given dictionary
            
            d dict: Key may be everything, Value is a number of states(cases) assotiated with the Key;
            
            Returns float: entropy value;
    '''
    a = np.array(d.values(), dtype=float)
    probs = a/np.sum(a)
    
    return -1*sum([(p*log(p)) for p in probs])
		

		
def get_simple_stat(array):
    return len(array), sum(array), array.mean(), np.median(array), min(array), max(array)


def counter2array(counter):
    return np.array(list(chain(*[[x[0]]*x[1] for x in counter.items()])))
	
		
def dict2simple_stat(counter):	
    return get_simple_stat(counter2array(counter))


def isinteger(s):
    '''Checks if provided string is convertible to integer'''
    try: 
        int(s)
        return True
    except ValueError:
        return False
            
def get_closest(a, n):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(a, n)
    if pos == 0:
        return a[0]
    if pos == len(a):
        return a[-1]
    before = a[pos - 1]
    after = a[pos]
    if after - n < n - before:
       return after
    else:
       return before
   
   
def hamming(s1, s2):
    """Calculate the Hamming distance between two strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def CDF(mylist, zerovalue=None):
    '''Creates CDF (cumulative distribution function) representation on basis of the provided iterable of floats/integers
    
        mylist list: element is float or tuple
        zerovalue int: if set this value will correspond to the CDF=0;
        
        Returns tuple: two elements: first is X-values of CDF, Second is Y-values of CDF 
    '''
    mylist.sort();
    norma = len(mylist)
    xvals, yvals = np.unique(mylist, return_index=True)
    if(zerovalue is None):
        yvals = np.array(list(yvals)[1:] + [norma])
    else:
        xvals = np.array([zerovalue] +list(xvals))
        yvals = np.array(list(yvals) + [norma])
    yvals = yvals/norma;
    return xvals, yvals


def lists2thresholds(l1, l2, reverse=False):
    '''For the two provided lists creates the dictionary, where the keys are all unique values in l1 and l2; and values are the numbers of l1 and l2 items which pass the threshold: item<=key (or item>=key if reverse==True);
    
        l1 list: elements are numbers
        l2 list: elements are numbers
        reverse bool: if set to True the logic changes to item>=key;
        
        Returns dictionary: keys are all unique values in l1 and l2; values 2-tuple: 1st is number of l1 items which pass the threshold: item<=key (or item>=key if reverse==True), 2nd is number of l1 items which pass the threshold: item<=key (or item>=key if reverse==True)
    '''
        
    l1.sort()
    l2.sort()
    x1, y1 = np.unique(l1, return_index=True)
    x2, y2 = np.unique(l2, return_index=True)
    norma1 = len(l1)
    norma2 = len(l2)
    y1 = np.append(y1[1:], norma1)
    y2 = np.append(y2[1:], norma2)
    d1 = dict(zip(x1, y1))
    d2 = dict(zip(x2, y2))
    uvals = np.append(x1, x2)
    uvals.sort()
    uvals = np.unique(uvals)
    
    res = {};
    cur1, cur2 = 0, 0;
    
    if(reverse):
        for threshold in uvals:
            res[threshold] = (norma1-cur1, norma2-cur2);
            cur1 = d1.get(threshold, cur1)
            cur2 = d2.get(threshold, cur2)
            
    else:
        for threshold in uvals:
            cur1 = d1.get(threshold, cur1)
            cur2 = d2.get(threshold, cur2)
            res[threshold] = (cur1, cur2);
        
    return res
        
        
        
def evaluate_increasing_trend(mylist, penalty):
    '''Evaluates how strong and stable is the increasing trend in the given list
        mylist list: list of numbers
        penalty float: if mylist[i+1] < mylist[i], the (mylist[i] - mylist[i+1])*penalty will be subtracted from the score. Must be bigger than 1.
        
        Returns float: evaluated score
    '''
    scores = [x[1]-x[0] for x in zip(mylist, mylist[1:])]
    return sum([x if x>0 else penalty*x for x in scores])


def evaluate_dereasing_trend(mylist, penalty):
    '''Evaluates how strong and stable is the increasing trend in the given list
        mylist list: list of numbers
        penalty float: if mylist[i+1] > mylist[i], the (mylist[i+1] - mylist[i])*penalty will be subtracted from the score. Must be bigger than 1.
        
        Returns float: evaluated score
    '''
    scores = [x[0]-x[1] for x in zip(mylist, mylist[1:])]
    return sum([x if x>0 else penalty*x for x in scores])


def find_best_trend(array2d, axis, penalty, increasing=True):
    if(increasing):
        evaluate = evaluate_increasing_trend
    else:
        evaluate = evaluate_dereasing_trend
    bestindex = np.argmax(np.apply_along_axis(evaluate, axis, array2d, penalty));
    if(axis == 0):
        return array2d[:,bestindex], bestindex
    else:
        return array2d[bestindex,:], bestindex
    
    
def find_elements_order(mylist):
    temp = sorted(mylist)     
    return [ (mylist.count(x)-1)/2 + temp.index(x) for x in mylist]
    
    
def get_accumulated_derivatives(xvalues, yvalues, window):
    """derivative averaged over a window of data"""
    return [(x[0] - x[1])/(x[2]-x[3]) for x in zip(yvalues[2*window:], yvalues, xvalues[2*window:], xvalues)]
    
        
def  smooth_with_averaging(vals, window):
    '''converts vals into the averages over the regions with flanks equal to window
            vals list: list of consecutive values
            window int: length of the flank. that is average is taken over window*2+1
    '''
    for w in range(window):
        yield np.mean(vals[:w*2+1])
    for el in sliding_window(vals, window*2+1):
        yield(np.mean(el))
    for w in range(window-1, -1, -1):
        yield np.mean(vals[-w*2-1:])
    
    
        
    
    
    
    
    
    
    


	
	
#testing section
if(__name__ == "__main__"):
    #mylist = [2,3,4, 1, 6, 1, 3, 7, 2]
    #print(mylist);
    #print(find_elements_order(mylist))
    
    a = list(range(20)) + [30]
    for el in smooth_with_averaging(a, 4):
        print(el)
