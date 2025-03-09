import os
import numpy as np
import math

import argparse

parser = argparse.ArgumentParser(description='xyz_2_gen')

# Arguments supported by the code.
parser.add_argument("--molecule1",      default='mol-1.xyz', help='mol-1.xyz')
parser.add_argument("--molecule2",      default='mol-2.xyz', help='mol-2.xyz')
parser.add_argument("--num_cell_input",       type=int,   default=3, help='3, 9, 10')
args        = parser.parse_args()
molecule1 = args.molecule1
molecule2 = args.molecule2
num_cell_input = args.num_cell_input


# compute distance using PBC, working only for orthorhombic cell?
def distance(x0, x1, cell):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * cell, delta - cell, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))



f1  = open(molecule1,      "r")

while True:
    tmp  = f1.readline()
    line = tmp.strip()
    if line == '': break

    tmp = tmp.split()
    natom1 = int(tmp[0])
    print ("number of atom ", natom1)
    
    tmp  = f1.readline()

    atomList1 = []
    xyz1 = np.zeros(shape=(natom1,3))
    for k in range(natom1):
        tmp = f1.readline().split()
        atomList1.append(tmp[0])
        xyz1[k,:] =  tmp[1:4]
    xyz1_center = np.mean(xyz1, axis=0)
f1.close()    

f2  = open(molecule2,      "r")

while True:
    tmp  = f2.readline()
    line = tmp.strip()
    if line == '': break

    tmp = tmp.split()
    natom2 = int(tmp[0])
    print ("number of atom ", natom2)

    tmp  = f2.readline().split()
    if num_cell_input==3:
        cell_abc = [float(x) for x in tmp]
        cell_abc = np.array(cell_abc)

    atomList2 = []
    xyz2 = np.zeros(shape=(natom2,3))
    for k in range(natom2):
        tmp = f2.readline().split()
        atomList2.append(tmp[0])
        xyz2[k,:] =  tmp[1:4]
    xyz2_center = np.mean(xyz2, axis=0)

f2.close()
print ("center of mass")
print (xyz1_center)
print (xyz2_center)

dist =  distance(xyz1_center, xyz2_center, cell_abc)
print ("dist=",dist)

f2 = open('TMP.xyz', "w")
f2.write("%-d\n" %(natom1+natom2))

if num_cell_input==3:
    for i in range(3):
        f2.write("  %15.9f" % cell_abc[i])
    f2.write("\n")

for i in range(natom1):
    f2.write("%-s %15.9f %15.9f %15.9f\n" %(atomList1[i], xyz1[i,0], xyz1[i,1], xyz1[i,2] ))
for i in range(natom2):
    f2.write("%-s %15.9f %15.9f %15.9f\n" %(atomList2[i], xyz2[i,0], xyz2[i,1], xyz2[i,2] ))

#f2.close()
#
