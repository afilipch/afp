'''library of classes and functions to extend and decorate functionality of numpy module'''

import copy;

import numpy as np



def merge_coordinates(c1,c2):
    '''Merges two area coordinates into one. For example two squares with coordinates [[1,5], [2,6]] and [[1,5], [4,9]] will be merged into one with coordinates [[1,5], [2,9]]
    
        c1 list: coordinates to merge. Elemenent is a list representing start and end of area at a certain dimension. 1st element of tuple is start(0-based inclusive), 2nd - end(exclusive)
        c2 list: coordinates to merge. Elemenent is a list representing start and end of area at a certain dimension. 1st element of tuple is start(0-based inclusive), 2nd - end(exclusive)
        
        Return list: elemenent is a list representing start and end of area at a certain dimension. 1st element of tuple is start(0-based inclusive), 2nd - end(exclusive)
    '''
    c = []
    for (s1,e1), (s2,e2) in zip(c1,c2):
        c.append([min(s1,s2), max(e1,e2)]);
    return c;	

    
def get_slice(a, coordinates):
    '''returns slice of numpy.array with given coordinates.
            
        a numpy.array: original array to produce slice from
        coordinates list: coordinates of slice. Elemenent is a list representing start and end of area at a certain dimension. 1st element of tuple is start(0-based inclusive), 2nd - end(exclusive)
        
        Return numpy.array: slice bounded by coordinates provided.
    '''
    expr = "a[" 
    for c in coordinates:
        expr += "%d:%d," % tuple(c);
    expr = expr[:-1] + "]";
    return eval(expr);
            
            
def one_side(coordinates, dimension, lookforward, start=True):
    '''extends coordinates in a given dimension
            
        coordinates list: coordinates to extend. Elemenent is a list representing start and end of area at a certain dimension. 1st element of tuple is start(0-based inclusive), 2nd - end(exclusive)
        dimension int: in this dimension coordinates are extended
        lookforward int: how far coordinates will are extended
        start bool: if True coordinates are extended in the left(bottom) direction. If False coordinates are extended in the right(upper) direction.
        
        Return list: extended coordinates. Elemenent is a list representing start and end of area at a certain dimension. 1st element of tuple is start(0-based inclusive), 2nd - end(exclusive)
    '''	
    nc = copy.deepcopy(coordinates);
    if(start==True):
        nc[dimension][1] = nc[dimension][0];
        nc[dimension][0] = nc[dimension][0] - lookforward;
    else:
        nc[dimension][0] = nc[dimension][1];
        nc[dimension][1] = nc[dimension][1] + lookforward;	
    return nc;
    
    
def keymax(a, key=lambda x: x):
    '''returns maximum(according to key function provided) elemenent of numpy.array
    
        a numpy.array: array to look for maximum element in
        key callable: function applied to the elements of array. The element with the maximum image(defined by the function) is returned.
        
        Return elemenent: maximum(according to key function provided) elemenent of numpy.array. The element with the maximum image(defined by the key function) is returned.
    '''	
    aiter = np.nditer(a)
    m = next(aiter);
    for c in aiter:
        if(key(c)>key(m)):
            m = c;
    return m;
    
    
def keymin(a, key=lambda x: x):
    '''returns the minimum(according to key function provided) elemenent of numpy.array
    
        a numpy.array: array to look for minimum element in
        key callable: function applied to the elements of array. The element with the minimum image(defined by the function) is returned.
        
        Return elemenent: the minimum(according to key function provided) elemenent of numpy.array. The element with the minimum image(defined by the key function) is returned.
    '''		
    aiter = np.nditer(a)
    m = next(aiter);
    for c in aiter:
        if(key(c)<key(m)):
            m = c;
    return m;	
    
    
def key_arg_max(a, key=lambda x: x):
    '''returns coordinates of the maximum(according to key function provided) elemenent of numpy.array
    
        a numpy.array: array to look for maximum element in
        key callable: function applied to the elements of array. Coordinates of the element with the maximum image(defined by the function) is returned.
        
        Return tuple: coordinates of the maximum(according to key function provided) elemenent of numpy.array. The element with the maximum image(defined by the key function) is returned.
    '''
    aiter = np.nditer(a);
    m = next(aiter);
    index = 0;
    for i,c in enumerate(aiter):
        if(key(c)>key(m)):
            m = c;
            index = i+1;
    return np.unravel_index(index, a.shape);
    
    
def key_arg_min(a, key=lambda x: x):
    '''returns coordinates of the minimum(according to key function provided) elemenent of numpy.array
    
        a numpy.array: array to look for minimum element in
        key callable: function applied to the elements of array. Coordinates of the element with the minimum image(defined by the function) is returned.
        
        Return tuple: coordinates of the minimum(according to key function provided) elemenent of numpy.array. The element with the minimum image(defined by the key function) is returned.
    '''	
    aiter = np.nditer(a);
    m = next(aiter);
    index = 0;
    for i,c in enumerate(aiter):
        if(key(c)<key(m)):
            m = c;
            index = i+1;
    return np.unravel_index(index, a.shape);

    
def counter2array(counter):
    pass;
    
#testing section	
if (__name__ == "__main__"):
    a = np.random.rand(3,5);
    #print a;
    #print 
    #print key_arg_min(a, lambda x: 1-(x-0.52)**2);


