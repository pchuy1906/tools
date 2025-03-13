import pandas as pd
import numpy as np

import argparse
parser = argparse.ArgumentParser(description='export FORCES, ENERGIES, and STRESS of file xyzfes')
# Arguments supported by the code.
parser.add_argument("--file_vDOS", default='power_spectrum_global.csv', help='file output from TRAVIS code')
parser.add_argument("--T", type=float, default=298.0, help='Temp. in K')
parser.add_argument("--nfu", type=int, default=20, help='number of atom per formula unit')
parser.add_argument("--mfu", type=float, default=12.0, help='formula unit mass (g/mol)')


args        = parser.parse_args()
file_vDOS           = args.file_vDOS
T                   = args.T
nfu                 = args.nfu
mfu                 = args.mfu


print ("-----------------------------------------------------------")
print ("Calculate the atomic energies")
print ("")
print ("read file_vDOS:", file_vDOS)

dataset = pd.read_csv(file_vDOS,delimiter=r';',skipinitialspace=True)
print (dataset.head())

names = ['#Wavenumber (cm^-1)', 'Spectrum (K*cm)']

freq      = dataset[names[0]].values
intensity = dataset[names[1]].values
freq = freq[1:]
intensity = intensity[1:]

inv_cm_to_K = 1.0/0.695

dw = freq[2]-freq[1]
nvib = 3.0 * float(nfu)
# normalize the vDOS by factor A, the new g(w) = A*intensity
A = nvib/( dw * sum(intensity))


# dimensionless frequency:
hbar_w__kB_T = freq * inv_cm_to_K/T
# dimensionless vDOS:
G = T/inv_cm_to_K * A * intensity

tmp = np.square(hbar_w__kB_T) * np.exp(hbar_w__kB_T) / np.square(np.exp(hbar_w__kB_T)-1.0) - 1.0
dx  = hbar_w__kB_T[1]-hbar_w__kB_T[0]
dCp = sum(tmp * G) * dx

kB = 1.380649 * np.power(10.0, -23)
nAvo = 6.02214076 * np.power(10.0, 23)
factor = kB/(mfu/nAvo)

print ("\nThe quantum correction for specific heat (J/g-K) is")
print (dCp*factor)
