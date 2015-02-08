h = 6.62606957 * 10 ** -34
#m = 1.67492735174 * 10 ** -27
m = 9.11 * 10**-31
L = 0.39 * 10**-9
convert  = 6.24150934 * 10**18

energylevels = []

def energy(n):
    return ((n*h/L)**2)/(8*m)*convert

for a in range(1,10000):
    energylevels.append(round(energy(a),4))


def energyassign(n):
    
