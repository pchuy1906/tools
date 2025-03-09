import os
import math
import numpy as np
from numpy import *

import argparse

parser = argparse.ArgumentParser(description='Calculate RMS')


# Arguments supported by the code.
parser.add_argument("--file_input", default='ENERGY.dat', help='file_input')


args        = parser.parse_args()

file_input   = args.file_input

print ("")
print ("Read input file ", file_input)
F = np.loadtxt(file_input)
print (F.shape, F.shape[0], F.shape[1] )


for i in range(F.shape[1]):
    Fave = sum(F[:,i])/float(F.shape[0])
    print Fave
    F[:,i] = F[:,i]-F[0,i]

np.savetxt('new_'+file_input, F, fmt='%15.9f %15.9f')

#RMS = 0
#ave_F_DFT = 0
#for i in range(F.shape[0]):
#    RMS += (F[i,0]-F[i,1])**2
#    ave_F_DFT += abs(F[i,0])
#RMS = np.sqrt ( RMS/ float(F.shape[0]) )
#ave_F_DFT = ave_F_DFT/float(F.shape[0])
#
#print "RMS=", RMS
#print "RMS_r=", RMS/ave_F_DFT

