from visual import *
from visual.graph import *
import numpy as np

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
maximum = 30

def oneDEnergy(n):
    energylevels = np.fromiter((((x*h/L)**2)/(8*m)*convert for x in range(1,n+1)),dtype=float)
    return energylevels
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

def boson(n,nthElvl,acc):
    '''
    n is the number of particles
    nthElvl is the nth energy level to go up to
    acc is integer of histogram(s)
    '''
    Elvl= oneDEnergy(nthElvl)
    #print Elvl
    econf= np.sort(configs(n, Elvl, particle="boson"))
    max_energy=Elvl[-1]+(n-1)*Elvl[0]
    #np.sort(econf)
    fcn1 = ghistogram(bins = np.linspace(econf[0],max_energy,acc))
    bound=True
    m=0
    while bound==True:
        if econf[m]>max_energy:
            bound=False
        else:
            m=m+1
            if m==econf.size:
                break
    fcn1.plot(data=econf[:m])
    return
    


boson(5,100,100)

