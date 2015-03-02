import numpy as np

def configs(x,elvl,particle="boson",out=None):
    """
    Generate configs for bosons or fermions.

    Parameters
    ----------
    num : positive integer
        Number of particles in the system
    Energylvls : 1-D array
        List of valid energy states of the system
    particle : "boson" or "fermion"
    out : 1-D array
        Array to put configs in

    Returns
    -------
    out : 1-D array
        array of total energys of all valid configurations

    """

    dtype = elvl.dtype  
    n=elvl.size #number of energy levels
    
    if particle=="boson":
        if out is None:
            out=None
            out = np.zeros(bosonconfigs(n,x), dtype=dtype)
        if x==1:
            for i in range(n):
                out[i]=elvl[i]
        if x>1:
            k=0 #index for end of last input added
            for m in range(n): #m is the energy level index
                end=k+bosonconfigs(n-m,x-1) #last energy level 
                configs(x-1, elvl[m:],particle,out[k:end])
                for i in range(k,end):
                    out[i]+= elvl[m]
                k=end
        return out

    if particle=="fermion":
        if out is None:
            out=None
            out = np.zeros(choose(n,x), dtype=dtype)
        if x==1:
            for i in range(n):
                out[i]=elvl[i]
        if x>1:
            k=0 #index for end of last input added
            for m in range(n): #m is the energy level index
                end=k+choose(n-(m+1),x-1) #last energy level 
                configs(x-1, elvl[m+1:],particle,out[k:end])
                for i in range(k,end):
                    out[i]+= elvl[m]
                k=end
        return out

    
def factorial(n):
    out=1
    for x in range(n):
        out=out*(x+1)
    return out

def bosonconfigs(numelvl,numpart):
    x=numpart
    n=numelvl
    out=choose(x+n-1,x)
    return out

def choose(choices,x):
    if choices>=x:
        out=int(factorial(choices)/((factorial(x)*factorial(choices-x))))
    else:
        out=0
    return out
