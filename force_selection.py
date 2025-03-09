import os
import numpy as np

fileforce = "COMPARE_F.dat"
force = np.loadtxt(fileforce)

dF=force[:,0]-force[:,1]
abs_dF = np.absolute(dF)

max_abs_dF = max(abs_dF)

f2 = open("force_selection.dat", "w")
if (max_abs_dF>1000.0):
    f2.write("%-d\n" %(0))
else:
    f2.write("%-d\n" %(1))
f2.close()

