# /usr/bin/python
'''library supporting logical rules generator based on False discovery rate'''

import sys
import os
import copy
from collections import defaultdict
from math import log;

import numpy as np

from afbio import numpy_extension
from afbio.numerictools import overlap_hyperrectangles;
from afbio.itertools_extension import  iterable_to_dict




class DimensionException(Exception):
	pass;	




class Grid(object):
    '''Grid represents multidimensional space containing analyzed data in a format of numpy.array. Each dimension corresponds to the particular attribute(weight, length, temperature). Each cell carries complex number. Real part is equal to a number of cases with attributes corresponding coordinates of the cell, Imaginary to the number of control ones. 
    
    Attributes:
        array numpy.array: the main attribute of the class. All data are stored here. Each dimension corresponds to the particular attribute(weight, length, temperature). Each cell carries complex number. Real part is equal to a number of cases with attributes corresponding coordinates of the cell, Imaginary to the number of control ones.
        
        encoding_table list: connects attributes values to the array coordinates. It is valuable for space compression. Element is dict (Key: attribute value, Value: coordinate). Inverse to the decoding_table
        
        decoding_table list: connects attributes values to the array coordinates. It is valuable for space compression. Element is dict (Key: coordinate, Value: attribute value). Inverse to the encoding_table
        
        attribute_names list: names of the dimensions of an array
    '''	
    def __init__(self, array, encoding_table, attribute_names=None):
        self.array = array
        self.encoding_table = encoding_table;
        self.shape = array.shape
        self.decoding_table = [];

        for d in encoding_table:
            self.decoding_table.append(dict([ (x[1], x[0]) for x in d.items()]))
                
        if(not attribute_names):
            self.attribute_names = ["x%d" % (i+1) for i in range(len(self.shape))]
        elif(len(attribute_names) == len(self.shape)):
            self.attribute_names = attribute_names;
        elif(len(attribute_names) > len(self.shape)):
            self.attribute_names = attribute_names[:len(self.shape)];
            sys.stderr.write('number of attribute names is more than number of dimensions, taken only first ones\n')
        else:
            self.attribute_names = attribute_names + ["x%d" % (i+1) for i in range(len(self.shape)-len(attribute_names))]
            sys.stderr.write('number of attribute names is less than number of dimensions, dummy names are assigned\n')
                    
            
    @classmethod		
    def from_dict(cls, signal, control, attribute_names=None):
        '''Constructs class instance from the given dictionaries
        
            signal dict: Based on real data, carries information about number of cases with particluar values of attributes. Key: tuple of attribute values, Value: number of cases corresponding to these attributes.
            control dict: Based on control data, carries information about number of cases with particluar values of attributes. Key: tuple of attribute values, Value: number of cases corresponding to these attributes.
            
            Return Grid: class instance
        '''
        #keys represent all existing combinations of attributes
        keys = list(signal.keys()) + list(control.keys());
        
        #check if dictionary key of signal and control elements are of the same length - which means that all elements can be placed in the one n-dimension grid
        ndims = set([len(x) for x in keys]);
        if(len(ndims) == 1):
            ndim = ndims.pop();
        else:
            raise DimensionException('all keys in provided dictionary have to be of equal length\n')
        
        #fill encoding_table
        encoding_table = list();
        for d in range(ndim):
            evs = sorted(set([x[d] for x in keys]));
            td = {};
            for cv, ev in enumerate(evs):
                td[ev] = cv;
            encoding_table.append(td);
        
        #fill array with signal(real) values
        array = np.zeros([len(x) for x in encoding_table], dtype=complex)
        for k, v in signal.items():
            key = []
            for d, kd in enumerate(k):
                key.append(encoding_table[d][kd])
            array[tuple(key)] += v
        
        #fill array with control(imaginary) values
        for k, v in control.items():
            key = []
            for d, kd in enumerate(k):
                    key.append(encoding_table[d][kd])
            array[tuple(key)] += complex(0, v)
                
        return cls(array, encoding_table, attribute_names=attribute_names);
            
            
    @classmethod		
    def from_iterable(cls, signal, control, entry='list', attributes=[], attribute_names=None):
        '''Constructs class instance from the given iterables
        
            signal iterable: Based on real data, carries all in a table-interpretable manner (rows: entries, columns: attributes)			
            control iterable: Based on control data, carries all in a table-interpretable manner (rows: entries, columns: attributes)
            
            entry list|object: type of iterable entry
                    if list: entry will be treated as an iterable with __getitem__ method. Attributes has to be list of integers or None
                    if object: entry will be treated as an object with __getattr__ method. Attributes has to be list of string or None
                    
            attributes list: attributes used for filtering. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)];
            attribute_names list: names of the attributes
            
            Return Grid: class instance
        '''
        sd = iterable_to_dict(signal, entry, 'count', attributes=attributes);
        cd = iterable_to_dict(control, entry, 'count', attributes=attributes);
        
        if(not attribute_names and entry=='object' and attributes):
            attribute_names = attributes;
        
        return cls.from_dict(sd, cd, attribute_names=attribute_names);
            
            
    @classmethod
    def from_csv(cls, csv, delimiter="\t"):
        '''Constructs class instance from the given csv table. Is used in testing mode. Method is able to construct only 2d space. Is supposed to read tables produced by 'to_csv' method
        
            csv string: path to csv table
            delimiter string: csv table delimiter
            
            Return Grid: class instance
        '''		
        def _converter(l):
            a = l.strip().split(delimiter)
            a = [x.replace("|", "+") + "j" for x in a]
            return delimiter.join(a);
                
        array = np.genfromtxt((_converter(x) for x in open(csv)), dtype=str, delimiter=delimiter)
        array = np.complex_(array)
        encoding_table = [dict([(x,x) for x in range(array.shape[0])]), dict([(x,x) for x in range(array.shape[1])])]
        return cls(array, encoding_table);
            
            
            
    def to_csv(self, csv="grid.csv", delimiter="\t"):
        '''Converts class instance into csv table. Is used in testing mode. Method is able to construct a table only from 2d space.
        
            csv string: path to csv table
            delimiter string: csv table delimiter
            
            Return Null
        '''			
        if(len(self.shape)!=2):
            raise DimensionException("Grid.to_csv() can be applied only for 2-dimension Grid objects\n")
        else:
            with open(csv,'w') as f:
                for row in self.array:
                    f.write(delimiter.join(["%d|%d" % (x.real, x.imag) for x in row])+"\n")
                                    
            
            
            
            
            
class Area(object):
    '''Respresent one square shaped area. Instances of this class are used in cluster extension proedure
    
    Attributes:
        coordinates list: coordinates of slice. Elemenent is a list representing start and end of area at a certain dimension. 1st element of tuple is start(0-based inclusive), 2nd - end(exclusive)
        
        array numpy.array: the main attribute of the class. All data are stored here. Each dimension corresponds to the particular attribute(weight, length, temperature). Each cell carries complex number. Real part is equal to a number of cases with attributes corresponding coordinates of the cell, Imaginary to the number of control ones.
        
        signal int: number of all items enclosed by the area
        control int: number of all control entries enclosed by the area
        fdr float: false discovery rate of items enclosed by area
            
    '''
    def __init__(self, coordinates, grid):
        self.coordinates = coordinates;
        self.array = numpy_extension.get_slice(grid.array, coordinates);
        cum = np.sum(self.array)
        self.signal = cum.real;
        self.control = cum.imag;
        if(self.signal):
            self.fdr = self.control/(self.control+self.signal);
        else:
            self.fdr = 1;
    

def cell_in_area(coordinates, area):
    '''Tests if the cell with given coordinates is inside the area
    
        coordinates tuple: coordinates of the cell of interest
        area Area: area(slice) of some numpy.array
        
        Return bool: True if cell with given coordinates is inside the area
    '''
    for c, (l,u) in zip(coordinates, area.coordinates):
        if(c<l or c>=u):
            return False;
    return True;
	
	
def intersect_areas(first_area, second_area, grid):
    '''Returns intersection of two areas
            
        first_area Area: an area(slice) of some numpy.array
        second_area Area: another area(slice) of some numpy.array
        grid Grid: Parental grid of the areas
        
        Returns Area: intersection of two areas. Returns None if there is no intersection
    '''
    coordinates = overlap_hyperrectangles(first_area.coordinates, second_area.coordinates)
    return Area(coordinates, grid);


	
	
class Cluster(object):
    '''Cluster object represents one logical conjunction("AND") rule. It contains area of the Grid corresponding to the rule(it means that upper and lowwer boundaries of the area can be translated to the rule s1<=x1<=e1 and s2<=x2<=e2 and etc)
    
    Attributes:
            grid Grid: Parental grid of the cluster area
            origin list: coordinates of the grid cell, which served as an origin of the cluster
            ID int: Unique id fo cluster. May correspond to the order cluster was generated
            extensions list: potential extensions to cluster area. Element is Area corresponding to the extension in one direction
            coordinates list: coordinates of slice. Elemenent is a list representing start and end of area at a certain dimension. 1st element of tuple is start(0-based inclusive), 2nd - end(exclusive)
            history list: history of extensions. Elemenent is extension coordinates(lisr of 2-element lists)
            nclusters list: list of negative clusters(negation rule). Represents rectangle areas inside the cluster with high false discovery rate
            support float: fraction of real cases inside the cluster
            rule string: logical rule representing cluster
            filter string: logical filter(ready to be passed into exec()) representing cluster
	'''
    def __init__(self, grid, origin, ID):
        self.grid = grid;
        self.origin = origin;
        self.ID = ID;
        
        self.coordinates = [];
        for p in origin:
                self.coordinates.append([p,p+1])
        self.history = [copy.deepcopy(self.coordinates)];
        
        self.extensions = [];
        self.nclusters = []
        self.support = 0;
        self.rule = ''
        self.filter_ = ''
            
		
    def __str__(self):
        coordinates = ''.join(["%s: %d<->%d\n" % (x[0], x[1][0], x[1][1]) for x in zip(self.grid.attribute_names, self.coordinates)]);
        return "cluster ID: %s\ncoordinates of origin: %s\nsignal of origin: %d\ncontrol of origin: %d\nsignal: %d\tcontrol: %d\nfdr: %1.4f\nsupport: %1.4f\ncluster coordinates:\n%s" % (self.ID, self.origin, self.grid.array[self.origin].real, self.grid.array[self.origin].imag, self.area.signal, self.area.control, self.area.fdr, self.support, coordinates)
		
		
    def _torule(self):
        single_rules = [];
        for d, (name, (start, end)) in enumerate(zip(self.grid.attribute_names, self.coordinates)):
            if(end >= self.grid.shape[d]):
                if(start):
                    sr = "%s>=%s" % (name, self.grid.decoding_table[d][start]);
                else:
                    continue;
            else:
                if(start):
                    sr = "%s<=%s<%s" % (self.grid.decoding_table[d][start], name, self.grid.decoding_table[d][end]);
                else:
                    sr = "%s<%s" % (name, self.grid.decoding_table[d][end]);
            single_rules.append(sr)
        self.rule = " and ".join(single_rules)
		
		
    def _tofilter_index(self, indices):
        single_rules = [];
        for d, (index, (start, end)) in enumerate(zip(indices, self.coordinates)):
            if(end >= self.grid.shape[d]):
                if(start):
                    sr = "x[%d]>=%s" % (index, self.grid.decoding_table[d][start]);
                else:
                    continue;
            else:
                if(start):
                    sr = "%s<=x[%d]<%s" % (self.grid.decoding_table[d][start], index, self.grid.decoding_table[d][end]);
                else:
                    sr = "x[%d]<%s" % (index, self.grid.decoding_table[d][end]);
            single_rules.append(sr)	
        self.filter_ = " and ".join(single_rules);


            
    def _tofilter_attribute(self, attributes):
        single_rules = [];
        for d, (attr, (start, end)) in enumerate(zip(attributes, self.coordinates)):
            if(end >= self.grid.shape[d]):
                if(start):
                    sr = "x.%s>=%s" % (attr, self.grid.decoding_table[d][start]);
                else:
                    continue;
            else:
                if(start):
                    sr = "%s<=x.%s<%s" % (self.grid.decoding_table[d][start], attr, self.grid.decoding_table[d][end]);
                else:
                    sr = "x.%s<%s" % (attr, self.grid.decoding_table[d][end]);
            single_rules.append(sr)
        self.filter_ = " and ".join(single_rules)
	
	
    def to_rule(self):
        self._torule();
        if(self.nclusters):
            nc_rules = []
            for nc in self.nclusters:
                nc._torule();
                nc_rules.append(nc.rule);
            self.rule = self.rule.join(["(", ") and (not ("]) + ")) and (not (".join(nc_rules)	 + "))"
        if(not self.rule):
            self.rule = 'spans the whole grid'
        return self.rule	
		
		
    def to_filter_index(self, indices):
        self._tofilter_index(indices);
        if(self.nclusters):
            nc_filters = []
            for nc in self.nclusters:
                nc._tofilter_index(indices);
                nc_filters.append(nc.filter_);
            self.filter_ = self.filter_.join(["(", ") and (not ("]) + ")) and (not (".join(nc_filters)	 + "))"
        if(not self.filter_):
            self.filter_ = 'True';
        return self.filter_
		
		
    def to_filter_attribute(self, attributes):
        self._tofilter_attribute(attributes);
        if(self.nclusters):
            nc_filters = []
            for nc in self.nclusters:
                nc._tofilter_attribute(attributes);
                nc_filters.append(nc.filter_);
            self.filter_ = self.filter_.join(["(", ") and (not ("]) + ")) and (not (".join(nc_filters)	 + "))"
        if(not self.filter_):
            self.filter_ = 'True'
        return self.filter_
            


    def get_extensions(self, lookforward, fdr_of_extension):
        '''Looks for possible extensions of the clusters area
        
                lookforward int: controls how far extension can go along any direction
                fdr_of_extension fdr: maximum false discovery rate allowed for extensions
                
                Adds new extensions to self.extensions
        '''
        for d, (s, e) in enumerate(self.coordinates):
            for lf in range(1, lookforward+1):
                if(s>=lf):
                    c = numpy_extension.one_side(self.coordinates, d, lf, start=True)
                    area = Area(c, self.grid);		
                    if(area.fdr<fdr_of_extension):
                        self.extensions.append(area);
                if(lf+e<=self.grid.shape[d]):
                    c = numpy_extension.one_side(self.coordinates, d, lf, start=False)
                    area = Area(c, self.grid);				
                    if(area.fdr<fdr_of_extension):
                        self.extensions.append(area);
        return True;

		
    def select_extensions(self, fit_function):
        '''Selects the best extension and extend cluster area with it
                
                fit_function: fitness function for extensions
                
                Extends self.coordinates in chosen direction. Set self.extensions ot an empty list. Updates self.history
        '''
        if(self.extensions):
            ext = max(self.extensions, key=fit_function);
            #print
            #print self.coordinates
            #print ext.signal, ext.control, ext.coordinates;	
            self.history.append(ext.coordinates);
            self.coordinates = numpy_extension.merge_coordinates(self.coordinates, ext.coordinates)
            self.extensions = [];
            return True;
        else:
            return False;
		
		
    def expand(self, fit_function, lookforward, fdr_of_extension):
        '''Expands cluster area until there are no possible extensions availible(with fdr <= fdr_of_extension)
        
                lookforward int: controls how far extension can go along any direction
                fdr_of_extension fdr: maximum false discovery rate allowed for extensions
                fit_function: fitness function for extensions
                
                Expands cluster area until there are no possible extensions availible(with fdr <= fdr_of_extension)
        '''
        expandable = True;
        while(expandable):
            self.get_extensions(lookforward, fdr_of_extension);
            expandable = self.select_extensions(fit_function)
        self.area = Area(self.coordinates, self.grid);
        #self.support = self.area.signal/(np.sum(grid.array).real+0.01)
		
		
    def get_nclusters(self, support, ncsupport, maxiter, fdr, lookforward, fit_function):
            
        for ncid in range(maxiter):
            ncorigin = numpy_extension.key_arg_max(self.area.array, key=fit_function);
            ncorigin = tuple([x[0]+x[1][0] for x in zip(ncorigin, self.coordinates)])
            #origin also has to meet fdr requirements
            if(self.grid.array[ncorigin].imag/max(self.grid.array[ncorigin].imag + self.grid.array[ncorigin].real, 1.0) < fdr):
                return False;
                    
            nc = Negative_cluster(self.grid, self, ncorigin, "%s:n%d" % (self.ID, ncid+1));
            nc.expand(fit_function, lookforward, fdr);
            
            nc.support = nc.area.control/self.area.control+0.01; 
            #we check for support of negative cluster both on level of parental cluster and parental grid
            if(nc.support > ncsupport and (nc.area.signal/self.area.signal+0.01)*self.support>support):
                sys.stderr.write("%s\n" % nc)
                self.nclusters.append(nc);
                nc.area.array[:] = 0;
            else:
                self.grid.array[ncorigin] = 0;
        return True;


		
class Negative_cluster(Cluster):
    def __init__(self, grid, cluster, origin, ID):
        super(Negative_cluster,self).__init__(grid, origin, ID)
        self.cluster = cluster;
        
    def get_extensions(self, lookforward, fdr_of_extension):
        '''Looks for possible extensions of the clusters area
        
            lookforward int: controls how far extension can go along any direction
            fdr_of_extension fdr: maximum false discovery rate allowed for extensions
            
            Adds new extensions to self.extensions
        '''
        for d, (s, e) in enumerate(self.coordinates):
            for lf in range(1, lookforward+1):
                if(s-lf>=self.cluster.coordinates[d][0]):
                    c = numpy_extension.one_side(self.coordinates, d, lf, start=True)
                    area = Area(c, self.grid);		
                    if(area.fdr>fdr_of_extension):
                        self.extensions.append(area);
                if(lf+e<=self.cluster.coordinates[d][1]):
                    c = numpy_extension.one_side(self.coordinates, d, lf, start=False)
                    area = Area(c, self.grid);
                    if(area.fdr>fdr_of_extension):
                        self.extensions.append(area);
        return True;
        

def get_rule(clusters):
    if(len(clusters)>1):
        return '(' + ') or ('.join(cl.to_rule() for cl in clusters) + ')'
    elif(clusters):
        return clusters[0].to_rule()
    else:
        return ""


def get_filter_index(clusters, indices):
    if(len(clusters)>1):
        return '(' + ') or ('.join(cl.to_filter_index(indices) for cl in clusters) + ')'
    else:
        return clusters[0].to_filter_index(indices)


def get_filter_attribute(clusters, attributes):
    if(len(clusters)>1):
        return '(' + ') or ('.join(cl.to_filter_index(attributes) for cl in clusters) + ')'
    else:
        return clusters[0].to_filter_index(attributes)
            
            
def apply_filter(signal, filter_):
    return eval('filter(lambda x: %s, %s)' % (filter_, signal))
            
            
def ff_fdr(area):
    '''fitness function for extensions'''
    return 1-area.fdr;
    
    
def ff_balanced(x):
    '''fitness function for extensions'''	
    if(type(x) == Area):
        return x.signal*(maxfdr-x.fdr);
    elif(type(x)==np.ndarray):
        if(x.real):
            return x.real*(maxfdr - x.imag/(x.imag+x.real))
        else:
            return maxfdr-1;
    
    
def nc_balanced(x):
    '''fitness function for extensions'''	
    if(type(x) == Area):
        return x.control*(x.fdr**2);
    elif(type(x)==np.ndarray):
        if(x.imag):
            return x.imag*((x.imag/(x.imag+x.real))**2)
        else:
            return 0;


def generate_clusters(grid, support = 0.01, maxiter = 100,  fdr=0.1, lookforward=10, fit_function=ff_balanced, ncsupport=0.1, nciter=0, ncfunction=nc_balanced):
    #maxfdr is used in fit_function, which can take only one argument(used as key function). So maxfdr has to be set as global
    global maxfdr
    maxfdr = fdr;
    
    #all generated clusters will be stored in 'clusters'
    clusters = [];
    total_signal = np.sum(grid.array).real+0.01;

    #we want to find an origin for a new cluster in that part of greed which is not occupied by any other cluster already generated
    free = Grid(np.copy(grid.array), grid.encoding_table, grid.attribute_names);

    for cid in range(maxiter):
        origin = numpy_extension.key_arg_max(free.array, key=fit_function);
        #origin also has to meet fdr requirement
        if(free.array[origin].imag/max(free.array[origin].imag + free.array[origin].real, 1.0) > fdr):
            return clusters
                
        cluster = Cluster(free, origin, "c%d" % (cid+1));
        cluster.expand(fit_function, lookforward, fdr);
        
        #test if cluster meet support requirements
        cluster.support = cluster.area.signal/total_signal
        if(cluster.support > support):
            sys.stderr.write("%s\n" % cluster);
            clusters.append(cluster);
            cluster.get_nclusters(support, ncsupport, nciter, fdr*2, lookforward, ncfunction)
            cluster.area.array[:] = 0;
            #cluster.to_filter_index([1,2])
        else:
            free.array[origin] = 0;

    return clusters;
    
    
    
    
def total_statistics(grid, clusters):
    signal = 0;
    for cl in clusters:
        signal += cl.area.signal - sum([x.area.signal for x in cl.nclusters]);
    control = 0;
    for cl in clusters:
        control += cl.area.control - sum([x.area.control for x in cl.nclusters]);
            
    support = signal/(np.sum(grid.array).real+0.01);
    
    if(signal):
        fdr = control/(control+signal);
    else:
        fdr = 1;
            
    return signal, control, support, fdr	


def lrg(signal, control, entry='list', attributes=[], attribute_names=[], support = 0.01, maxiter = 10,  fdr=0.1, lookforward=10, fit_function=ff_balanced, ncsupport=0.1, nciter=0, ncfunction=nc_balanced):

    grid = Grid.from_iterable(signal, control, entry=entry, attributes=attributes, attribute_names=attribute_names);
    clusters = generate_clusters(grid, support = support, maxiter = maxiter,  fdr=fdr, lookforward=lookforward, fit_function=fit_function, ncsupport=ncsupport, nciter=nciter, ncfunction=ncfunction);
    
    if(not clusters):
        return None, None, None
    
    rule = get_rule(clusters)
    if(entry == 'list'):
        if(not attributes):
            attributes= list(range(len(attribute_names)))
        lrg_filter = get_filter_index(clusters, attributes)
    else:
        lrg_filter = get_filter_attribute(clusters, attributes)
            
    signal_total, control_total, support_total, fdr_total = total_statistics(grid, clusters);
    
    log_message = "\nfilter applied: %s\n\nrule generated: %s\n\nnumber of instances passed filter: %d\nnumber of control instances passed filter: %d\nfraction of instances passed filter: %1.5f\nestimated FDR: %1.5f\n" % (lrg_filter, rule, signal_total, control_total, support_total, fdr_total)
            
    return lrg_filter, rule, log_message
		
	

#def reassign_grid(clusters, grid):
    #for cluster in clusters:
            #cluster.grid = grid;
            #cluster.area = Area(cluster.coordinates, grid);



#testing section
if(__name__ == "__main__"):
    from random import randint;
    signal = [(randint(1,6), randint(1,6), randint(1,6), randint(2,8)) for _ in range(10000)]
    signal += [(randint(12,16), randint(12,16), randint(14,20), randint(10,14)) for _ in range(10000)]
    control = [(randint(1,20), randint(1,20), randint(1,20), randint(1,20)) for _ in range(2000)]
    
    
    filter_, rule = lrg(signal, control, entry='list', attributes=[], attribute_names=None, support = 0.01, maxiter = 10,  fdr=0.1, lookforward=10, fit_function=ff_balanced, ncsupport=0.1, nciter=0, ncfunction=nc_balanced)
    filtered = apply_filter(signal, filter_)
    
    #print filter_
    #print
    #print len(filtered)
    
	
	
	
	
	

