import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='Normalize to range [0,maxval]')

# Arguments supported by the code.
parser.add_argument("--file_input", default='data.dat', help='file_input')
parser.add_argument("--maxval", type=float, default=1.0, help='max value (min value = 0)')

args        = parser.parse_args()

file_input   = args.file_input
maxval   = args.maxval

print ("")
print ("Read input file ", file_input)
F = np.loadtxt(file_input)

nline = F.shape[0]
Fmin = min(F[:,1])
Fmax = max(F[:,1])

print ("the number of line: ", nline)
print ("Fmin, Fmax=", Fmin, Fmax )

for i in range(nline):
    F[i,1] = (F[i,1]-Fmin)/(Fmax-Fmin)*maxval

file_output = "normalized_" + file_input
np.savetxt(file_output, F, fmt='%15.9f %15.9f')

