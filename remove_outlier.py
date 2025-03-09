import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='remove outlier')


# Arguments supported by the code.
parser.add_argument("--file_input", default='COMPARE.dat', help='file_input')


args        = parser.parse_args()

file_input   = args.file_input

print ("")
print ("Read input file ", file_input)
F = np.loadtxt(file_input, usecols=[0,1])
print ("The size of the array:")
print (F.shape, F.shape[0])

print ()
print ("after removing outliers")
newF = F[np.abs(F[:,1]-F[:,0]) < 200]
print ("The size of the array:")
print (newF.shape, newF.shape[0])

np.savetxt("data_final.dat", newF)

