import numpy as np
from visual import *
from visual.graph import *
#functions
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
        for j in range(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

#variables
energylevels=np.fromiter((x for x in range(10)),dtype=int)


def bosonconfigs(n,Earray,Emax, out=None):
    """
    Generate the configurations possible bosons 

    Parameters
    ----------
    n: integer
        number of bosons
    Earray: ndarray
        ordered array of possible energy states of a single boson
    Emax: float
        maximum energy of the configurations

    Returns
    -------
    out: ndarray
        1-D array of shape m containing the energy of M possible configurations
        of bosons with total energy less than Emax.

    Example
    -------
    >>> bosonconfigs(3,np.fromiter(x for x in range(5)),3.5)
    array([ 0,
            1,
            2,
            3,
            1,
            2,
            3,
            1,
            2,
            3,
            3,
            2,
            3,
            3])

    Note this array:
    array([[0, 0, 0],
           [0, 0, 1],
           [0, 0, 2],
           [0, 0, 3],
           [0, 1, 0],
           [0, 1, 1],
           [0, 1, 2],
           [0, 2, 0],
           [0, 2, 1],
           [1, 0, 0],
           [1, 0, 1],
           [1, 0, 2],
           [1, 1, 0],
           [1, 1, 1],
           [1, 2, 0],
           [2, 0, 0],
           [2, 0, 1],
           [2, 1, 0]])
    """
    foo=Emax-Earray[0]*(n-1)
    x=0
    bar=Earray[0]
    while bar<foo and Earray[x]:
       x=x+1
       bar=Earray[x]
       if out is None:
           m=(x+1)**n
           out=np.zeros([m],dtype=float)

    if n>1:
       t=0 
    if n==1:
        for i in range(x):
            out[x]=out[x]+Earray[x]
    
    return out

def lessindex(array, maximum):
    index=0
    amax=array[0]
    while amax<maximum:
       index=index+1
       if index==array.size:
           break
       amax=array[index]
    return index-1
            
def unoptomizedbosonconfigs(n,Earray,Emax,res):
    Earraymax=Emax-Earray[0]*(n-1)
    x=lessindex(Earray,Earraymax)
    configs=cartesian([Earray[0:x] for y in range(n)])
    Econfigs=np.ndarray(configs.size/n,dtype=float)
    for i in range(configs.size/n):
        for j in configs[i]:
            Econfigs[i]+=j
    a = (max(Econfigs) - min(Econfigs))/res
    fnc1 = ghistogram(bins = arange(0,Emax,a), color = color.red)
    fnc1.plot(data=([x for x in Econfigs if x<Emax]))

h = 6.62606957 * 10 ** -34      #Plank's constant
#m = 1.67492735174 * 10 ** -27   #this is mass of neutron
m = 9.11 * 10**-31             #this is mass of electron
L = 0.39 * 10**-9               #size of box
convert  = 6.24150934 * 10**18
maximum = 30


energylevels = np.fromiter((((x*h/L)**2)/(8*m)*convert for x in range(maximum+1)),dtype=float)
#this creates the entire table of energy levels as a single list
#unoptomizedbosonconfigs(5,energylevels,1000,300)

threedconfig=cartesian([energylevels for x in range(3)])
threedenergylevel=np.fromiter((np.sum(threedconfig[x]) for x in range(threedconfig.size/3)),dtype=float)
threedenergylevel.sort()
a = (max(threedenergylevel) - min(threedenergylevel))/1000
fnc2 = ghistogram(bins = arange(0,max(threedenergylevel)), color = color.red)
fnc2.plot(data=threedenergylevel)
