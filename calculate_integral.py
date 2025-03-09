import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='Calculate integral')


# Arguments supported by the code.
parser.add_argument("--file_input", default='gr_H_H.dat', help='file_input')

args        = parser.parse_args()

file_input   = args.file_input

print ("")
print ("Read input file ", file_input)
F = np.loadtxt(file_input)
print (F.shape, F.shape[0])

fx = []
for i in range(F.shape[0]):
    fx.append(F[i,1]* F[i,0]*F[i,0])

print ("Estimated integral:")
print ( (F[1,0]-F[0,0])* sum(fx) )

n = F.shape[0]
print (F[n/2:n/2+10])

f2 = open('gr_r_r.dat', "w")
for i in range(F.shape[0]):
    f2.write("%15.9f %15.9f\n" %( F[i,0], fx[i] ))
f2.close
