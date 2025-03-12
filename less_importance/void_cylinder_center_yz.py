import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import random

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--Rpore", type=float, default=10, help='pore radius (in Angstrom)')
parser.add_argument("--natom_per_molecule", type=int, default=1, help='number of atom per molecule')


args               = parser.parse_args()
file_xyz           = args.file_xyz
Rpore              = args.Rpore
natom_per_molecule = args.natom_per_molecule

print ("")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

cell_xyz = np.zeros(shape=(3,3))
tmp = f.readline()
tmp = tmp.split()
try:
    # format: cell_xyz[1,:] cell_xyz[2,:] cell_xyz[3,:]
    k = 0
    for k1 in range(3):
        for k2 in range(3):
            cell_xyz[k1,k2] = float(tmp[k])
            k += 1
except:
    try:
        # format: NON_ORTHOR cell_xyz[1,:] cell_xyz[2,:] cell_xyz[3,:]
        k = 1
        for k1 in range(3):
            for k2 in range(3):
                cell_xyz[k1,k2] = float(tmp[k])
                k += 1
    except:
        cell_xyz = np.zeros(shape=(3,3))
        # Orthor cell a, b, c
        for k1 in range(3):
            cell_xyz[k1,k1] = float(tmp[k1])

print (cell_xyz)
atomList = []
xyz = np.zeros(shape=(natom,3))
for k in range(0,natom):
    tmp = f.readline()
    tmp = tmp.split()
    atomList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
f.close

asym = np.array(atomList)
asym_unique = np.unique(asym)
asym_unique = np.sort(asym_unique)
asym_list = ' '.join(asym_unique)

f2 = open( "file_out.xyz", "w")
f2.write("%s\n" %( "@natom@"))
#f2.write("%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %d %s\n" %( cell_xyz[0,0], cell_xyz[0,1], cell_xyz[0,2], cell_xyz[1,0], cell_xyz[1,1], cell_xyz[1,2], cell_xyz[2,0], cell_xyz[2,1], cell_xyz[2,2], len(asym_unique), asym_list ))
f2.write("%15.9f %15.9f %15.9f %d %s\n" %( cell_xyz[0,0], cell_xyz[1,1], cell_xyz[2,2], len(asym_unique), asym_list ))


nmolecule = natom/natom_per_molecule

def check_remove(xyz, Rpore, cell_xyz):
    center = (cell_xyz[0,:] + cell_xyz[1,:] + cell_xyz[2,:]) * 0.5
    dxyz  = center - xyz
    dxyz  = np.array(dxyz)
    R2 = dxyz[1]*dxyz[1] + dxyz[2]*dxyz[2] 
    if (R2 < Rpore*Rpore):
        return 1
    else:
        return 0

for i in range(nmolecule):
    iremove = 0
    for j in range(natom_per_molecule):
        k = i*natom_per_molecule + j
        iremove += check_remove(xyz[k], Rpore, cell_xyz)
    if (iremove <= natom_per_molecule/2):
        for j in range(natom_per_molecule):
            k = i*natom_per_molecule + j
            f2.write("%4s %15.9f %15.9f %15.9f\n" %( atomList[k], xyz[k,0], xyz[k,1], xyz[k,2] ))
f2.close()


