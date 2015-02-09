#this is test code to ensure that github is working properly

from numpy import histogram
from math import factorial
from itertools import combinations
from visual import *
from visual.graph import *

h = 6.62606957 * 10 ** -34
#m = 1.67492735174 * 10 ** -27
m = 9.11 * 10**-31
L = 0.39 * 10**-9
convert  = 6.24150934 * 10**18
maximum = 100

energylevels = []

def energy(n):
    return ((n*h/L)**2)/(8*m)*convert

for a in range(1,maximum):
    energylevels.append(round(energy(a),4))


def fermion(number):
    energycount = []
    #generate some list of all possible combinations of energy levels (imcomplete)
    fermionlist = list(combinations(energylevels,int(number)))
    #calculate total energy of each configuration
    for config in fermionlist:
        total = 0
        for i in config:
            total = total + energy(i)
        energycount.append(total)
    #return number of configurations in energy range
    a = (max(energycount) - min(energycount))/300
    fnc1 = ghistogram(bins = arange(min(energycount), max(energycount), int(a)), color = color.red)
    fnc1.plot(data=energycount)
    return histogram(energycount,bins = 300,weights = None, density = False)
    
        
