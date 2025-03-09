import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import argparse
parser = argparse.ArgumentParser(description='group same molecules to one xyz file')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')

args        = parser.parse_args()
file_xyz           = args.file_xyz

print ("-----------------------------------------------------------")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"rt")

def number_CHNO(atom_order, molecule_order, molecule_natom):
    res = np.array([0,0,0,0])
    for id_tmp in range(len(molecule_order)):
        atom_loc = np.where(atom_order == molecule_order[id_tmp])
        res[atom_loc] = molecule_natom[id_tmp]
    return res

def iname(molecule_order, molecule_natom):
    nlen = len(molecule_order)
    molecule_name = ""
    for i in range(nlen):
        aaa = molecule_order[i] + str(molecule_natom[i]) + "_"
        molecule_name = molecule_name + aaa
    return molecule_name


istruc = 0
atom_order = np.array(['C','H','N','O'])

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break
    #v, e = [float(x) for x in line.split()[:2]]

    natom = int(tmp)
    
    box = np.zeros(shape=(3,3))

    tmp = f.readline().split()
    for k in range(3):
        box[k,k] = float(tmp[k])
    energy = tmp[3]

    atomList = []
    xyz = np.zeros(shape=(natom,3))
    fxyz = np.zeros(shape=(natom,3))

    for k in range(natom):
        tmp = f.readline().split()
        atomList.append(tmp[0])
        xyz[k,0], xyz[k,1], xyz[k,2] = float(tmp[1]),float(tmp[2]),float(tmp[3])
        fxyz[k,0],fxyz[k,1],fxyz[k,2] = float(tmp[4]),float(tmp[5]),float(tmp[6])
    molecule_order, molecule_natom = np.unique(atomList, return_counts=True)
    molecule_name = iname(molecule_order, molecule_natom)
    print (molecule_name)

    istruc += 1
    fname = molecule_name + ".xyz"
    f2  = open(fname ,"a+")
    f2.write("%1d\n" %( natom ))
    f2.write("%15.9f %15.9f %15.9f %s\n" %( box[0,0], box[1,1], box[2,2], energy))
    for k in range(natom):
        f2.write("%s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n" %(atomList[k], xyz[k,0], xyz[k,1], xyz[k,2], fxyz[k,0],fxyz[k,1],fxyz[k,2]))


f.close

