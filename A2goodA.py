import os
import numpy as np
from math import floor

import argparse
parser = argparse.ArgumentParser(description='xyz_2_lammps for SHOCK')
# Arguments supported by the code.
parser.add_argument("--file_xyz",                         default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--num_cell_input",       type=int,   default=3, help='3, 9, 10')

args        = parser.parse_args()
file_xyz        = args.file_xyz
num_cell_input  = args.num_cell_input


f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

cell_3_3 = np.zeros(shape=(3,3))
tmp = f.readline().split()

if num_cell_input==3:
    for k in range(3):
        cell_3_3[k,k] = float(tmp[k])
if num_cell_input==9:
    cell_3_3[0,0] = float(tmp[0])
    cell_3_3[0,1] = float(tmp[1])
    cell_3_3[0,2] = float(tmp[2])

    cell_3_3[1,0] = float(tmp[3])
    cell_3_3[1,1] = float(tmp[4])
    cell_3_3[1,2] = float(tmp[5])

    cell_3_3[2,0] = float(tmp[6])
    cell_3_3[2,1] = float(tmp[7])
    cell_3_3[2,2] = float(tmp[8])

if num_cell_input==10:
    cell_3_3[0,0] = float(tmp[1])
    cell_3_3[1,1] = float(tmp[5])
    cell_3_3[2,2] = float(tmp[9])

myList = []
xyz = np.zeros(shape=(natom,3))
for k in range(natom):
    tmp = f.readline()
    tmp = tmp.split()
    myList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
f.close

numrows = xyz.shape[0]
numcols = xyz.shape[1]
print ("size of xyz:", numrows, numcols)

print ("")
print "Angstrom to Crystal calculation..."
Cxyz = np.dot(xyz, np.linalg.inv(cell_3_3))

for k1 in range(0,numrows):
    for k2 in range(0,numcols):
        tmp = Cxyz[k1][k2]-floor(Cxyz[k1][k2])
        Cxyz[k1][k2] = tmp

print "Crystal to Angstrom calculation..."
Axyz = np.dot(Cxyz, cell_3_3)

print ("")
print "Output: <Axyz.xyz>"
f2 = open('Axyz.xyz', "w")

f2.write("%-d\n" %(natom))

print (cell_3_3)
for ixyz in range(3):
    f2.write("  %15.9f" % cell_3_3[ixyz,ixyz])
f2.write("\n")

for k in range(0,natom):
    f2.write("%4s %15.9f %15.9f %15.9f\n" %( myList[k], Axyz[k,0], Axyz[k,1], Axyz[k,2] ))
f2.close
