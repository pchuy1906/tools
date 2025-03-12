cm_inv_2_kcal_mol = 0.0028591

import numpy as np
file_VDOS = "power_spectrum_global.csv"

print ("Read file VDOS (power_spectrum_global.csv)")
data = np.genfromtxt(file_VDOS, delimiter=';')

freq = data[:,0]
inte = data[:,1]
print
print ("Calculate the integral")
ZPE = sum(freq * inte) / sum(inte)
ZPE *= cm_inv_2_kcal_mol

print ("The zero point energy is", ZPE, " kcal/mol")

