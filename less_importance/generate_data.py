import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='generate the data from fluctuate data')

# Arguments supported by the code.
parser.add_argument("--file_input1", default='data1.dat', help='file_input')
parser.add_argument("--file_input2", default='data2.dat', help='file_input')
parser.add_argument("--file_output", default='data3.dat', help='file_output')

args        = parser.parse_args()

file_input1   = args.file_input1
file_input2   = args.file_input2
file_output   = args.file_output

print ("")
print ("Read input file1 ", file_input1)
F1 = np.loadtxt(file_input1)
nline = F1.shape[0]
nhalf = nline/2
mean1 = np.mean(F1[nhalf:,1])
print ("Read input file2 ", file_input2)
F2 = np.loadtxt(file_input2)
F2 = F2[1:,:]
mean2 = np.mean(F2[:,1])

print ("mean1 = ", mean1)
print ("mean2 = ", mean2)

F2[:,0] = F2[:,0] + F1[-1,0]
F2[:,1] = F2[:,1] + mean1 - mean2

F = np.vstack([F1, F2])

F[F<0] = 0

np.savetxt(file_output, F, fmt='%15.9f %15.9f')

