# /usr/bin/python
'''library supporting advanced handling, processing, conversion of iterables'''

import sys
import os
import itertools
from collections import defaultdict

import numpy

class LengthException(Exception):
    pass;
    
    
class ArgumentException(Exception):
    pass;	

    
    
    
def pop_first(iterable):
    '''pops first element of any iterable, regardless does it has __getitem__ method or not'''
    if(hasattr(iterable, 'next')):	
        return iterable.next();
    else:
        return iterable.pop(0);

            
def iterable_of_objects_to_counts_dict(iterable, attributes=[]):
    '''Returns dictionary build on iterable. Key: tuple of certain values of attributes, Value: number of objects in iterable with values of attributes equal to the ones in Key.
    
        iterable iterable: any iterable of objects		
        attributes list: attributes used to pull objects. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)]
        if None all attributes will be used
        
        Returns dictionary: Key: tuple of certain values of attributes, Value: number of objects in iterable with values of attributes equal to the ones in Key.
    '''
    d = defaultdict(int);	
    
    try:
        if(hasattr(iterable, 'next')):	
            first = iterable.next();
            if(not attributes):
                attributes += vars(first).keys()			
            k = tuple([getattr(first, attr) for attr in attributes])
            d[k]+=1;	
        else:
            first = iterable[0];
            if(not attributes):
                attributes += vars(first).keys()
    except(StopIteration):
        return d
    
    for el in iterable:
        k = tuple([getattr(el, attr) for attr in attributes])
        d[k]+=1;
    return d;
    
    
def iterable_of_lists_to_counts_dict(iterable, indices=[]):
    '''Returns dictionary build on iterable. Key: tuple of certain values of attributes(corresponding to indices), Value: number of objects in iterable with values of attributes equal to the ones in Key.

    iterable iterable: any iterable of objects		
    indices list: indices used to pull objects. For example if only index 3 is provided all entries with the entry[3] will be gathered in one dictionary entry with Key=[(entry[3])]
    if None all attributes will be used

    Returns dictionary: Key: tuple of certain values of attributes(corresponding to indices), Value: number of objects in iterable with values of attributes equal to the ones in Key.
    '''
    d = defaultdict(int);
    if(indices):
        for el in iterable:	
            k = tuple([el[i] for i in indices])
            d[k]+=1;
    else:
        for el in iterable:
            k = tuple([x for x in el])
            d[k]+=1;
    return d;
    
	
def iterable_of_objects_to_list_dict(iterable, attributes=[]):
    '''Returns dictionary build on iterable. Key: tuple of certain values of attributes, Value: list of objects in iterable with values of attributes equal to the ones in Key.
    
        iterable iterable: any iterable of objects		
        attributes list: attributes used to pull objects. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)]
        if None all attributes will be used
        
        Returns dictionary: Key: tuple of certain values of attributes, Value: list of objects in iterable with values of attributes equal to the ones in Key.
    '''
    d = defaultdict(list);

    try:
        if(hasattr(iterable, 'next')):	
            first = iterable.next();
            if(not attributes):
                attributes += vars(first).keys()			
            k = tuple([getattr(first, attr) for attr in attributes])
            d[k].append(first);	
        else:
            first = iterable[0];
            if(not attributes):
                attributes += vars(first).keys()
    except(StopIteration):
        return d
                    
    for el in iterable:
        k = tuple([getattr(el, attr) for attr in attributes])
        d[k].append(el);
    return d;
	
	
def iterable_of_lists_to_list_dict(iterable, indices=[]):
    '''Returns dictionary build on iterable. Key: tuple of certain values of attributes(corresponding to indices), Value: list of objects in iterable with values of attributes equal to the ones in Key.
    
        iterable iterable: any iterable of objects		
        indices list: indices used to pull objects. For example if only index 3 is provided all entries with the entry[3] will be gathered in one dictionary entry with Key=[(entry[3])]
        if None all attributes will be used
        
        Returns dictionary: Key: tuple of certain values of attributes(corresponding to indices), Value: list of objects in iterable with values of attributes equal to the ones in Key.
    '''	
    d = defaultdict(list);
    
    try:
        if(hasattr(iterable, 'next')):	
            first = iterable.next();
            length = len(first);		
            if(not indices):
                indices += range(length)			
            k = tuple([first[i] for i in indices])
            d[k].append(first);
        else:
            first = iterable[0];
            length = len(first);
            if(not indices):
                indices += range(length)
    except(StopIteration):
        return d
            
    for el in iterable:
        if(len(el)!=length):
            raise LengthException('elements in iterable are of different length\n')
        else:
            k = tuple([el[i] for i in indices])
            d[k].append(el);
    return d;
	

def iterable_to_dict(iterable, entry, mode, attributes=[]):
    '''Returns dictionary build on iterable. Key: tuple of certain values of attributes, Value: number of objects in iterable with values of attributes equal to the ones in Key.
    
        iterable iterable: any iterable of objects	
        
        entry list|object: type of iterable entry
            if list: entry will be treated as an iterable with __getitem__ method. Attributes has to be list of integers or None
            if object: entry will be treated as an object with __getattr__ method. Attributes has to be list of string or None
        
        mode list|count: 
            if list: entries with the same attributes will pulled in a list
            if count: number of entries with the same attributes will be stored
                
        attributes list: attributes used to pull objects. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)]
        OR indices used to pull objects. For example if only index 3 is provided all entries with the entry[3] will be gathered in one dictionary entry with Key=[(entry[3])]
        if None all attributes will be used
        
        Returns dictionary: Key: tuple of certain values of attributes, Value: number or list of objects in iterable with values of attributes equal to the ones in Key.
    '''
    if(entry == 'object' and mode == 'count'):
        return iterable_of_objects_to_counts_dict(iterable, attributes=attributes);
    elif(entry == 'object' and mode == 'list'):
        return iterable_of_objects_to_list_dict(iterable, attributes=attributes);
    elif(entry == 'list' and mode == 'count'):
        return iterable_of_lists_to_counts_dict(iterable, indices=attributes);
    elif(entry == 'list' and mode == 'list'):
        return iterable_of_lists_to_list_dict(iterable, indices=attributes);
    else:
        raise ArgumentException('Wrong arguments were passed to \'iterable_to_dict\'\nentry has to set to \'list\' or \'object\'\nmode has to set to \'list\' or \'count\'\n')
		
	

	
	
	
	
def cmp_indices(l1, l2, indices):
    for i in indices:
        r = cmp(l1[i], l2[i]);
        if(r):
            return r;
    else:
        return 0;
            
            
def cmp_attributes(obj1, obj2, attributes):
    for attr in attributes:
        r = cmp(getattr(obj1, attr), getattr(obj2, attr));
        if(r):
            return r;
    else:
        return 0;
            
                    
def flatten(iterable):
    return list(itertools.chain.from_iterable(iterable))
    
    
def counter2list(counter):
    '''Flatten collections.Counter-like dictionary into  \'original list\', that is: Counter(counter2list(counter)) == counter
    
            counter dict: collections.Counter-like dictionary to flatten
            
            Returns list: list of counter keys, duplicated \'corresponding\' value times
    '''	
            
    l = [];
    for k, v in counter.iteritems():
        l += [k]*v
    return l;
    
    
def median(lst):
    return numpy.median(numpy.array(lst))


def rolling_window(lst, window_size):
    '''Generates consecutive windows of size [window_size] from lst iterable'''
    it = iter(lst)
    win = [it.next() for cnt in xrange(window_size)] # First window
    yield win
    for e in it: # Subsequent windows
        win[:-1] = win[1:]
        win[-1] = e
        yield win


def smooth_values(lst, window_size):
    '''Generates sums of consecutive widows of size [window_size] from lst iterable'''
    it = iter(lst)
    win = [it.next() for cnt in xrange(window_size)] # First window
    aggregate = sum(win)
    for _ in win:
        yield aggregate
            
    for e in it: # Subsequent windows
        aggregate = aggregate - win[0] + e;
        win[:-1] = win[1:]
        win[-1] = e
        yield aggregate
            
            
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))
    

#testing section
if(__name__ == "__main__"):
    pass;
	
	
	
	