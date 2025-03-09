import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands


import argparse
parser = argparse.ArgumentParser(description=' select data from 1d array')
# Arguments supported by the code.
parser.add_argument("--file_input", default='file.dat', help='file input 1d array')
parser.add_argument("--distance", type=float, default=5.0, help='distance to select item')


args         = parser.parse_args()
file_input   = args.file_input
distance     = float(args.distance)


print ("")
print ("read file input:", file_input)
print ("")


F12 = np.loadtxt(file_input)
arr_1d = F12[:,0]

res = []

curr_posi = min(arr_1d)
res.append((F12[0,0],F12[0,1]))
for k in range(1,len(arr_1d)):
    if (abs(curr_posi-arr_1d[k]) > distance):
	res.append((F12[k,0],F12[k,1])) 
	curr_posi = F12[k,0]

print (len(res))

np.savetxt("new_"+file_input, res)


#np.savetxt('boxXYZ', boxXYZ, fmt='%15.9f %15.9f %15.9f')

