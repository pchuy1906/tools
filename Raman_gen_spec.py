import glob
raman_files = glob.glob('vasp_raman_*.dat')

import numpy as np
npoint = 800
freq = np.linspace(100.0, 1000.0, num=npoint)
Inte = np.zeros(shape=(npoint))

T = 300
# Kelvin to cm-1 
K_2_inv_cm = 0.695

# hbar unit eV/cm-1
hbar = 0.00012398425731484318

# kB unit eV/K
kB = 8.617333262 * 0.00001

def adding_peak(T, freq, Inte, file):
    with open(file, 'r') as f2:
	last_line = f2.readlines()[-1]
    tmp = last_line.split()
    if (tmp[0] != "#"):
        f = float(tmp[1])
        I = float(tmp[4])

        beta = 1.0/(kB*T)
        f_quantum = beta*hbar*f/(1.0-np.exp(-beta*hbar*f))
        f_quantum = 1.0

        FWHM = T * K_2_inv_cm
        FWHM = 10
        FWHM = FWHM*0.5
        sig = FWHM/2.355

        Inte += I*np.exp(-(freq-f)**2/(2.0* sig* sig))* f_quantum

for file in raman_files:
    adding_peak(T, freq, Inte, file)

combined = np.vstack((freq, Inte)).T
np.savetxt('RAMAN.dat', combined, fmt='%15.9f %15.9f')

from scipy.signal import find_peaks
id_peaks, _ = find_peaks(Inte, height=0)
#xplot = freq[id_peaks]
#yplot = Inte[id_peaks]
combined = np.vstack((freq[id_peaks], Inte[id_peaks])).T
np.savetxt('list_freq.dat', combined, fmt='%15.9f %15.9f')

import matplotlib.pyplot as plt
plt.plot(freq, Inte)
#plt.plot(xplot, yplot,"x")

plt.show()
