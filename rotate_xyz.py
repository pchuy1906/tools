import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import argparse
parser = argparse.ArgumentParser(description='rotate file XYZ')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--iatom1", type=int, default=1, help='atom-1')
parser.add_argument("--iatom2", type=int, default=2, help='atom-2')
parser.add_argument("--iatom3", type=int, default=3, help='atom-3')
parser.add_argument("--rot_ang", type=float, default=90.0, help='rotate angle')

args        = parser.parse_args()
file_xyz    = args.file_xyz
iatom1      = args.iatom1
iatom2      = args.iatom2
iatom3      = args.iatom3
rot_ang     = args.rot_ang

print ("")
print ("Doing rotation for file:")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

print ("")
natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)
tmp = f.readline()

myList = []
xyz = np.zeros(shape=(natom,3))
for k in range(0,natom):
    tmp = f.readline()
    tmp = tmp.split()
    myList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
f.close

print ("")
print ("Atom 1 and atom 2 are the axis to rotate. Details are:")
print ("Atom 1: ", " ID = ", iatom1, myList[iatom1-1], xyz[iatom1-1,:])
print ("Atom 2: ", " ID = ", iatom2, myList[iatom2-1], xyz[iatom2-1,:])
print ("Atom 3 is selected for the plane of the molecule. Details are:")
print ("Atom 3: ", " ID = ", iatom3, myList[iatom3-1], xyz[iatom3-1,:])
print ("")


print ("building the axis:")
x = xyz[iatom2-1,:] - xyz[iatom1-1,:]
x = x/np.linalg.norm(x)
yp= xyz[iatom3-1,:] - xyz[iatom1-1,:]
z = np.cross(x,yp)
z = z/np.linalg.norm(z)
y = np.cross(z,x)

print ("")
print ("The threee vectors are:")
print (x)
print (y)
print (z)
print ("")

box = np.vstack((x, y, z))
print ("merge to a matrix")
print (box)
print ("")

print ("Calculating the internal coordinates:")
Cxyz = np.dot(xyz, np.linalg.inv(box))

print ("")
print ("The angle to be rotated is:", rot_ang)
print ("")
print ("Calculating the new matrix")


ang_rad = rot_ang/180.0*3.14159265359

xp =  x
yp =  y*cos(ang_rad) + z*sin(ang_rad)
zp = -y*sin(ang_rad) + z*cos(ang_rad)
newbox = np.vstack((xp, yp, zp))

print (newbox)

print ("")
print ("Calculating the new coordinates")
newxyz = np.dot(Cxyz, newbox)
#print (type(newxyz), newxyz.shape)
tran = -newxyz[iatom1-1,:] + xyz[iatom1-1,:]
for k in range(newxyz.shape[0]):
    newxyz[k,:] = newxyz[k,:] + tran[:]

print (newxyz)

f2 = open("rotated.xyz", "w")
f2.write("%1d\n" %( newxyz.shape[0] ))
f2.write("%1s\n" %( "" ))

f3 = open("xyz.rotated", "w")
for k in range(newxyz.shape[0]):
    f2.write("%1s %15.9f %15.9f %15.9f\n" %(myList[k],newxyz[k,0],newxyz[k,1],newxyz[k,2] ))
    f3.write("%1s %15.9f %15.9f %15.9f\n" %(myList[k],newxyz[k,0],newxyz[k,1],newxyz[k,2] ))
f2.close()
f3.close()


