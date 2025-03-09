import os
import numpy as np
from numpy import *
import commands

import argparse

parser = argparse.ArgumentParser(description='Extract info to XYZFES file')

# Arguments supported by the code.
parser.add_argument("--file_qeoutput", default='output', help='file_qeoutput')
parser.add_argument("--ioption", type=int, default=2, help='1-vcrelax   2-relax/scf')

args       = parser.parse_args()
file_qeoutput   = args.file_qeoutput
ioption         = args.ioption

# unit conversion
bohr2angstrom = 0.52917720859
Ry2Ha = 0.500
Ry2kcalmol = 313.75470835207074
kbar2GPa = 0.100

# initialize variables:
celldm = 0.0

box = np.zeros(shape=(3,3))

filename = "DFT.xyzfes"
f2 = open(filename, "w")

with open(file_qeoutput, 'rb') as f:
    while True:
        line = f.readline()
        if not line:
            break

        # read the number of atoms
        if "number of atoms/cell" in line:
           natom = line.split()[4]
           natom = int(natom)

        # read the cell-parameters
        if "celldm(1)" in line:
           celldm = line.split()[1]
           celldm = float(celldm)
        if "a(1)" in line:
           tmp = line.split()[3:6]
           for i in range(3):
              box[0,i] = float(tmp[i]) * celldm * bohr2angstrom
        if "a(2)" in line:
           tmp = line.split()[3:6]
           for i in range(3):
              box[1,i] = float(tmp[i]) * celldm * bohr2angstrom
        if "a(3)" in line:
           tmp = line.split()[3:6]
           for i in range(3):
              box[2,i] = float(tmp[i]) * celldm * bohr2angstrom

        # read the atomic coordinates
        # case-1
        if "Cartesian axes" in line:
           line = f.readline()
           line = f.readline()
           atomList = []
           xyz = np.zeros(shape=(natom,3))
           for i in range(natom):
              line = f.readline()
              tmp = line.split()
              atomList.append(tmp[1])
              xyz[i,0] =  float(tmp[6]) * celldm * bohr2angstrom
              xyz[i,1] =  float(tmp[7]) * celldm * bohr2angstrom
              xyz[i,2] =  float(tmp[8]) * celldm * bohr2angstrom
        # case-2
        if "ATOMIC_POSITIONS (angstrom)" in line:
           atomList = []
           xyz = np.zeros(shape=(natom,3))
           for i in range(natom):
              line = f.readline()
              tmp = line.split()
              atomList.append(tmp[0])
              xyz[i,0] =  float(tmp[1])
              xyz[i,1] =  float(tmp[2])
              xyz[i,2] =  float(tmp[3])

        # read the energy
        if "!    total energy" in line:
           tmp = line.split()
           energy =  float(tmp[4]) * Ry2kcalmol

        # read the forces
        if "Forces acting on atoms" in line:
           line = f.readline()
           fxyz = np.zeros(shape=(natom,3))
           for i in range(natom):
              line = f.readline()
              tmp = line.split()
              fxyz[i,0] =  float(tmp[6]) * Ry2Ha
              fxyz[i,1] =  float(tmp[7]) * Ry2Ha
              fxyz[i,2] =  float(tmp[8]) * Ry2Ha

        # read the stress tensor
        if "total   stress  (Ry/bohr**3)" in line:
           line = f.readline()
           tmp = line.split()
           sxx = float(tmp[3]) * kbar2GPa
           sxy = float(tmp[4]) * kbar2GPa
           sxz = float(tmp[5]) * kbar2GPa
           line = f.readline()
           tmp = line.split()
           #syx = float(tmp[3]) * kbar2GPa
           syy = float(tmp[4]) * kbar2GPa
           syz = float(tmp[5]) * kbar2GPa
           line = f.readline()
           tmp = line.split()
           #szx = float(tmp[3]) * kbar2GPa
           #szy = float(tmp[4]) * kbar2GPa
           szz = float(tmp[5]) * kbar2GPa
           
           # write the xyz file
           f2.write("%6d\n"%(natom))
           f2.write("%12s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n"%("NON_ORTHO",   box[0,0],box[0,1],box[0,2],   box[1,0],box[1,1],box[1,2],   box[2,0],box[2,1],box[2,2],  sxx,syy,szz,sxy,sxz,syz, energy))
           for k in range(0,natom):
               f2.write("%6s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n" %( atomList[k], xyz[k,0], xyz[k,1], xyz[k,2], fxyz[k,0], fxyz[k,1], fxyz[k,2] ))
f2.close() 
