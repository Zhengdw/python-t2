import numpy as np
from math import *
from visual import *
from visual.graph import *


def energy2(n):
    return ((n*h/L)**2)/(8*m)*convert

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

def configs(x,elvl,particle="boson",out=None):
    """
    Generate configs for bosons or fermions.

    Parameters
    ----------
    x : positive integer
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

h = 6.62606957 * 10 ** -34      #Plank's constant
#m = 1.67492735174 * 10 ** -27   #this is mass of neutron
m = 9.11 * 10**-31             #this is mass of electron
L = 0.39 * 10**-9               #size of box
convert  = 6.24150934 * 10**18
maximum = 20
kb = 1.3806488 * 10**-23

energylevels = np.fromiter((((x*h/L)**2)/(8*m)*convert for x in range(maximum+1)),dtype=float)
#this creates the entire table of energy levels as a single list

def ThreeDEnergy(n):
    out=np.empty(n**3)
    index=0
    for i in range(1,n+1):
        for j in range(1,n+1):
            for k in range(1,n+1):
                out[index]=((i*h/L)**2+(j*h/L)**2+(k*h/L)**2)/(8*m)*convert
                index=index+1
    return np.sort(out)

def fermion(n):
    energycount = []
    energy = ThreeDEnergy(maximum)
    fermionlist = configs(n,energy,'fermion')
    for config in np.nditer(fermionlist):
        if config < energy2(maximum):
            energycount.append(config)
    #return number of configurations in energy range
    a = (max(energycount) - min(energycount))/200
    fnc1 = ghistogram(bins = arange(min(energycount), max(energycount), int(round(a))), color = color.red)
    fnc1.plot(data=energycount)
    hist, binedges = np.histogram(energycount,bins = 100,weights = None, density = False)
    return (hist,binedges)
   
    
def boltzfit(xvalues,yvalues,degree):
    xlist = []
    ylist = []
    for i in range(len(xvalues)-1):
        xlist.append((xvalues[i] + xvalues[i+1])/(2*convert)) #Average energy in joules
    for j in range(len(yvalues)):
        ylist.append(kb*log(yvalues[j])) #Convert to boltzmann entropy
    return np.polyfit(xlist,ylist,degree)

def gibbsfit(xvalues,yvalues,degree):
    xlist = []
    for i in range(len(xvalues)-1):
        xlist.append((xvalues[i] + xvalues[i+1])/(2*convert)) #Average energy in joules
    for j in range(len(yvalues)): #generate gibbs configurations
        if j > 0:
            yvalues[j] = yvalues[j-1] + yvalues[j]
    for j in range(len(yvalues)):
        ylist.append(kb*log(yvalues[j])) #Convert to entropy
    return np.polyfit(xlist,ylist,degree)

def entropy(n,degree,particle ='boson',method ='gibbs'):
    if particle == 'boson':
        data = boson(n)
    if particle == 'fermion':
        data = fermion(n)
    if method == 'boltzmann':
        return (boltzfit(data[1],data[0],degree))
    if method == 'gibbs':
        return (gibbsfit(data[1],data[0],degree))        
                   


    
    



