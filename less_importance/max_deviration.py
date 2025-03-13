import os
import math
import numpy as np
import pandas as pd
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='calculate maximum deviation')

# Arguments supported by the code.
parser.add_argument("--file_input", default='external_pressure.dat', help='file_input')

args        = parser.parse_args()

file_input   = args.file_input

print ("")
print ("Read input file ", file_input)
F = np.loadtxt(file_input)

Ndata = F.shape[0]
print ("shape of the file:")
print (F.shape)

dF = F[:,0]-F[:,1]
dF = np.abs(dF)

print ("min and max differences:")
print (np.min(dF), np.max(dF))

AA = np.where(dF>100)
print (AA)
print (F[AA,:])

np.savetxt('ID_to_compare.dat', AA, fmt='%i', delimiter="\n")

