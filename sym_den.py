import os
import math
import numpy as np

import argparse

parser = argparse.ArgumentParser(description='calculate average')


# Arguments supported by the code.
parser.add_argument("--file_input", default='density.dat', help='file_input')
args      = parser.parse_args()
file_input   = args.file_input

F = np.loadtxt(file_input)
X = F[:,0]
Y = F[:,1]
X = np.array(X)
Y = np.array(Y)

Ndata = F.shape[0]
print ("shape of the file:")
print (F.shape)

#res = np.vstack((ave, std)).T
#newfile = "results_ave_std.dat"
#np.savetxt(newfile, res, fmt='%15.9f %15.9f')

indices = np.where(Y > 0)[0]
idmin = indices[0]
idmax = indices[-1]
idmin = idmin - 2
idmax = idmax + 2
print (idmin, idmax)

idmiddle_real = 0.5 * float(idmin+idmax)
print (idmiddle_real)
idmiddle_int = int(idmiddle_real)

for i in range(idmin, idmiddle_int+1):
    i1 = i
    i2 = idmax - i + idmin
    F[i1,0] = 0.5*( X[i1]-X[i2])
    F[i2,0] = 0.5*(-X[i1]+X[i2])
    F[i1,1] = 0.5*(Y[i1]+Y[i2])
    F[i2,1] = 0.5*(Y[i1]+Y[i2])
    print (i, i1, i2)

newfile = "sym_" + file_input
idx=0
np.savetxt(newfile, F[idmin-idx:idmax+idx,:], fmt='%15.9f %15.9f')

