import os
import numpy as np
from numpy import *
import commands

import argparse

parser = argparse.ArgumentParser(description='stretch atom in XYZ file')

# Arguments supported by the code.
parser.add_argument("--iatom1", type=int, default=1, help='index of atom as origin')
parser.add_argument("--iatom2", type=int, default=2, help='index of atom as stretching')
parser.add_argument("--file_input", default='file.xyz', help='file_xyz')
parser.add_argument("--delta", type=float, default=0.2, help='to be stretch by a step delta Angstrom')
parser.add_argument("--rmax", type=float, default=1.95, help='to be stretch until rmax Angstrom')


args       = parser.parse_args()
file_input   = args.file_input
iatom1       = args.iatom1
iatom2       = args.iatom2
delta        = args.delta
rmax         = args.rmax


def makeXYZ(filename, xyz, myList):
    natom = xyz.shape[0]
    f2 = open(filename, "w")
    f2.write("%6d\n"%(natom))
    f2.write("%6s\n"%("comment"))
    for k in range(0,natom):
        f2.write("%6s %15.9f %15.9f %15.9f\n" %( myList[k], xyz[k,0], xyz[k,1], xyz[k,2] ))
    f2.close()


f  = open(file_input ,"r")
# read line 1: natom
tmp = f.readline()
natom = int(tmp)
# read line 2: comment
tmp = f.readline()
#tmp = tmp.split()
#tmp = np.array(tmp)
#boxXYZ = [float(x) for x in tmp]
#boxXYZ = np.array(boxXYZ)
#boxXYZ = boxXYZ.reshape((3, 3))
# read xyz
myList = []
xyz = np.zeros(shape=(natom,3))
for k in range(0,natom):
    tmp = f.readline()
    tmp = tmp.split()
    myList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])

# make XYZ files
r12 = xyz[iatom2-1,:]-xyz[iatom1-1,:]
d12 = sqrt(r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2])

#print (r12, d12)
istruc = 0

filename =  str(istruc) +".xyz"
makeXYZ(filename, xyz, myList)

while (d12+delta < rmax):
   istruc += 1
   xyz[iatom2-1,:] = xyz[iatom1-1,:] + r12 * (d12 + delta)/d12
   filename =  str(istruc) +".xyz"
   makeXYZ(filename, xyz, myList)

   r12 = xyz[iatom2-1,:]-xyz[iatom1-1,:]
   d12 = sqrt(r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2])
   #print (r12, d12)


