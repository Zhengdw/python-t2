import numpy as np
from math import factorial
from itertools import combinations
from visual import *
from visual.graph import *

h = 6.62606957 * 10 ** -34
#m = 1.67492735174 * 10 ** -27   #this is mass of neutron
m = 9.11 * 10**-31             #this is mass of electron
L = 0.39 * 10**-9
convert  = 6.24150934 * 10**18
maximum = 100

energylevels = np.fromiter((((x*h/L)**2)/(8*m)*convert for x in range(100)),dtype=float)
#this creates the entire table of energy levels as a single list

def energy(n):
    return ((n*h/L)**2)/(8*m)*convert

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
    a = (max(energycount) - min(energycount))/300
    fnc1 = ghistogram(bins = arange(min(energycount), max(energycount), int(round(a))), color = color.red)
    fnc1.plot(data=energycount)
    return np.histogram(energycount,bins = 300,weights = None, density = False)

# BTW, Jonathan, this isn't working for nuetron masses. Something to do with rounding
