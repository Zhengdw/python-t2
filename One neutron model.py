import numpy as np
from math import factorial
from itertools import combinations
from visual import *
from visual.graph import *

h = 6.62606957 * 10 ** -34      #Plank's constant
#m = 1.67492735174 * 10 ** -27   #this is mass of neutron
m = 9.11 * 10**-31             #this is mass of electron
L = 0.39 * 10**-9               #size of box
convert  = 6.24150934 * 10**18
maximum = 100

energylevels = np.fromiter((((x*h/L)**2)/(8*m)*convert for x in range(1,maximum)),dtype=float)
#this creates the entire table of energy levels as a single list

def energy(n):
    return ((n*h/L)**2)/(8*m)*convert

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype
 
    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out
'''
#def fermion(n):
    config = []
    energycount = []
    for i in range(n):
        config.append(i+1)
    while sum(config) < 3*max(energylevels) - 3:
        total = 0
        for part in config:
            total = total + part
        energycount.append(total)
        cont = True
        while cont:
            if config[n-1] != maximum:
                config[n-1] = config[n-1] + 1
            else:
                
'''       
        #I need to find some way to have the configuration update properly
        #Is there some way to make a recursive if loop?
    

def boson(number):
    energyconfigs=cartesian([energylevels[0:20] for x in range(number)])
    print (energyconfigs[1])
    


def fermion(number):
    energycount = []
    #generate some list of all possible combinations of energy levels (incomplete)
    fermionlist = list(combinations(energylevels,int(number)))
    #calculate total energy of each configuration
    for config in fermionlist:
        total = 0
        for i in config:
            total = total + energy(i)
        energycount.append(total)
    #return number of configurations in energy range
    maxi = 10**7
    a = (max(energycount) - min(energycount))/200
    fnc1 = ghistogram(bins = arange(min(energycount), max(energycount), int(round(a))), color = color.red)
    fnc1.plot(data=energycount)
    #return np.histogram(energycount,bins = 100,weights = None, density = False)




