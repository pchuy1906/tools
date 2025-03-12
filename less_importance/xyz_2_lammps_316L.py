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
parser.add_argument("--ishock", type=int, default=0, help='0-no shock, 1-shock')
parser.add_argument("--iSTIB", type=int, default=0, help='0-no STIBs, 1-STIBs')
parser.add_argument("--dSTIB", type=float, default=20.0, help='the wide of the STIBs')

args        = parser.parse_args()
file_xyz    = args.file_xyz
ishock      = args.ishock
iSTIB       = args.iSTIB
dSTIB       = args.dSTIB

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
for k in range(3):
    box[k,k] = float(tmp[k])

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
print (syms, counts_syms, ntype)

dz = 0.0
tz = 0.0

if (ishock == 1):
    deps = 0.01
    zmin = min(xyz[:,2])
    zmax = max(xyz[:,2])
    print ("zmin, zmax = ", zmin, zmax)
    tz = zmin-deps
    dz = zmax-zmin+deps-box[2,2]+ deps/2.0


f2 = open("new_data.lammps", "w")
f2.write("%1s\n" %( "# Position data file" ))
f2.write("%1s\n" %( "" ))
f2.write("%1d %1s\n" %( natom, "atoms" ))
f2.write("%1d %1s\n" %( ntype, "atom types" ))
f2.write("%1s\n" %( "" ))
f2.write("%15.9f %15.9f %1s\n" %( 0.0, box[0,0], "xlo xhi" ))
f2.write("%15.9f %15.9f %1s\n" %( 0.0, box[1,1], "ylo yhi" ))
Lz = box[2,2]+dz
f2.write("%15.9f %15.9f %1s\n" %( 0.0, Lz, "zlo zhi" ))
f2.write("%1s\n" %( "" ))
f2.write("%1s\n" %( "Masses" ))
f2.write("%1s\n" %( "" ))

#data_sym = ["C","H","O","N"]
#data_mass= [12.0107, 1.0079, 15.999, 14.0064]

#data_sym = ["N","H"]
#data_mass= [14.0064, 1.0079]

data_sym = ["Fe","Ni","Cr"]
data_mass= [55.845, 58.6934, 51.9961]

#data_sym = ["Fe"]
#data_mass= [55.845]


for k in range(len(data_sym)):
    print (k+1, data_mass[k])
    f2.write("%1d %15.9f\n" %( k+1, data_mass[k] ))

f2.write("%1s\n" %( "" ))
f2.write("%1s\n" %( "Atoms" ))
f2.write("%1s\n" %( "" ))

STIB1 = "group  stib1 id "
STIB2 = "group  stib2 id "

for k in range(0,natom):
    for isym in range(len(data_sym)):
        if (myList[k] == data_sym[isym]):
            atype = isym + 1
    #f2.write("%1d %1d %1d %15.9f %15.9f %15.9f\n" %( k+1, atype, 0, xyz[k,0], xyz[k,1], xyz[k,2]-1.0*tz ))
    f2.write("%1d %1d %15.9f %15.9f %15.9f\n" %( k+1, atype, xyz[k,0], xyz[k,1], xyz[k,2]-1.0*tz ))
    if (xyz[k,2]-1.0*tz < dSTIB): STIB1 += " " + str(k+1)
    if (xyz[k,2]-1.0*tz > Lz-dSTIB): STIB2 += " " + str(k+1)

f2.close()

if (iSTIB==1):
    outF = open("STIB1.txt", "w")
    outF.write(STIB1)
    outF.close()
    
    outF = open("STIB2.txt", "w")
    outF.write(STIB2)
    outF.close()

