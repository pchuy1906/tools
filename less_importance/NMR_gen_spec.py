file="sigma.dat"

import numpy as np

chemical_shifts = np.loadtxt(file)

minx = min(chemical_shifts)
maxx = max(chemical_shifts)
print ("chemical shift ranges:", minx, maxx)

npoint = 500
deps = 5.0
freq = np.linspace(minx-deps, maxx+deps, num=npoint)
Inte = np.zeros(shape=(npoint))

FWHM = 0.5
sig = FWHM/2.355

for k in range(len(chemical_shifts)):
    f = chemical_shifts[k]
    Inte += 1.0*np.exp(-(freq-f)**2/(2.0* sig* sig))

combined = np.vstack((freq, Inte)).T
np.savetxt('spectrum.dat', combined, fmt='%15.9f %15.9f')

#import matplotlib.pyplot as plt
#plt.plot(freq, Inte)
#plt.plot(xplot, yplot,"x")
#plt.show()
