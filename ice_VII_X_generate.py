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
parser.add_argument("--OH_dist", type=float, default=1.0, help='OH bond')

args        = parser.parse_args()
file_xyz     = args.file_xyz
OH_dist      = args.OH_dist


print ("read fileXYZ:", file_xyz)
f2 = open("file_new.xyz","w")
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)
f2.write("%d" %(natom))

tmp = f.readline()
f2.write("\n%s" %(tmp))

myList = []
xyz = np.zeros(shape=(natom,3))
for k in range(natom):
    tmp = f.readline()
    tmp = tmp.split()
    myList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
    if (k%3==2):
        f2.write("%4s %15.9f %15.9f %15.9f\n" %( myList[k], xyz[k,0], xyz[k,1], xyz[k,2] ))
        OH1 = xyz[k-1,:]-xyz[k,:]
        rOH1= np.sqrt(np.dot(OH1,OH1))
        OH2 = xyz[k-2,:]-xyz[k,:]
        rOH2= np.sqrt(np.dot(OH2,OH2))
        nH1 = xyz[k,:] + OH1 * OH_dist/rOH1
        nH2 = xyz[k,:] + OH2 * OH_dist/rOH2
        f2.write("%4s %15.9f %15.9f %15.9f\n" %( myList[k-1], nH1[0], nH1[1], nH1[2] ))
        f2.write("%4s %15.9f %15.9f %15.9f\n" %( myList[k-2], nH2[0], nH2[1], nH2[2] ))


f.close()
f2.close()



