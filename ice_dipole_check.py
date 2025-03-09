import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')

args        = parser.parse_args()
file_xyz     = args.file_xyz

print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

tmp = f.readline()

myList = []
xyz = np.zeros(shape=(natom,3))

dipX = 0.0
dipY = 0.0
dipZ = 0.0

for k in range(natom):
    tmp = f.readline()
    tmp = tmp.split()
    myList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
    if (k%3==2):
        dipX += -2.0*xyz[k,0]
        dipY += -2.0*xyz[k,1]
        dipZ += -2.0*xyz[k,2]
    else:
        dipX += 1.0*xyz[k,0]
        dipY += 1.0*xyz[k,1]
        dipZ += 1.0*xyz[k,2]
f.close

print ("Final dipole:%5.2f %5.2f %5.2f" % (dipX, dipY, dipZ) )


