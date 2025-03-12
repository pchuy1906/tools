import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import argparse
parser = argparse.ArgumentParser(description='xyz_2_lammps for SHOCK')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')

args        = parser.parse_args()
file_xyz    = args.file_xyz

print ("")
print ("Orthohombic cell case")
print ("")
print ("need to check variable data_sym and data_mass")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

box = np.zeros(shape=(3,3))
tmp = f.readline()
tmp = tmp.split()
box[0,:] = [float(x) for x in tmp[0:3]]
box[1,:] = [float(x) for x in tmp[3:6]]
box[2,:] = [float(x) for x in tmp[6:9]]


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

syms, counts_syms = np.unique(myList, return_counts=True)
ntype = len(syms)

f2 = open("data.lammps", "w")
f2.write("%1s\n" %( "# Position data file" ))
f2.write("%1s\n" %( "" ))
f2.write("%1d %1s\n" %( natom, "atoms" ))
f2.write("%1d %1s\n" %( ntype, "atom types" ))
f2.write("%1s\n" %( "" ))
f2.write("%15.9f %15.9f %1s\n" %( 0.0, box[0,0], "xlo xhi" ))
f2.write("%15.9f %15.9f %1s\n" %( 0.0, box[1,1], "ylo yhi" ))
f2.write("%15.9f %15.9f %1s\n" %( 0.0, box[2,2], "zlo zhi" ))
f2.write("%1s\n" %( "" ))
f2.write("%1s\n" %( "Masses" ))
f2.write("%1s\n" %( "" ))

data_sym = ["Si","O"]
data_mass= [28.0855, 15.999]

for k in range(len(data_sym)):
    print (k+1, data_mass[k])
    f2.write("%1d %15.9f\n" %( k+1, data_mass[k] ))

f2.write("%1s\n" %( "" ))
f2.write("%1s\n" %( "Atoms" ))
f2.write("%1s\n" %( "" ))

for k in range(0,natom):
    for isym in range(len(data_sym)):
        if (myList[k] == data_sym[isym]):
            atype = isym + 1
    #f2.write("%1d %1d %1d %15.9f %15.9f %15.9f\n" %( k+1, atype, 0, xyz[k,0], xyz[k,1], xyz[k,2] ))
    f2.write("%1d %1d %15.9f %15.9f %15.9f\n" %( k+1, atype, xyz[k,0], xyz[k,1], xyz[k,2] ))
f2.close()
