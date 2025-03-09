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
parser.add_argument("--OH_bond", type=float, default=1.0, help='OH_bond')

args        = parser.parse_args()
file_xyz     = args.file_xyz
OH_bond      = args.OH_bond

def Cell_XYZ_ABC(cellXYZ):
    a1= cellXYZ[0,:]
    a2= cellXYZ[1,:]
    a3= cellXYZ[2,:]
    a = math.sqrt(dot(a1, a1))
    b = math.sqrt(dot(a2, a2))
    c = math.sqrt(dot(a3, a3))
    alp = math.acos(dot(a2, a3)/(b*c))*180.0/pi
    bet = math.acos(dot(a1, a3)/(a*c))*180.0/pi
    gam = math.acos(dot(a1, a2)/(a*b))*180.0/pi
    return np.array([a,b,c,alp,bet,gam])

def Cell_ABC_XYZ(cellABC):
    tmp = np.zeros(shape=(3,3))
    tmp[0,0] = cellABC[0]
    tmp[1,0] = cellABC[1]*cos(cellABC[5])
    tmp[1,1] = cellABC[1]*sin(cellABC[5])
    tmp[2,0] = cellABC[2]*cos(cellABC[4])
    tmp[2,1] = cellABC[2]*cos(cellABC[3])*sin(cellABC[5])-((cellABC[2]*cos(cellABC[4])-cellABC[2]*cos(cellABC[3])*cos(cellABC[5]))/tan(cellABC[5]))
    tmp[2,2] = sqrt((cellABC[2])**2 -(tmp[2,0])**2 - (tmp[2,1])**2)
    return tmp

print ("")
print ("cell_9")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

natom = f.readline()
natom = int(natom)
print ("the number of atom:", natom)

box = np.zeros(shape=(3,3))
tmp = f.readline()
tmp = tmp.split()
tmp = np.array(tmp)

box[0,:] = [float(x) for x in tmp[0:3]]
box[1,:] = [float(x) for x in tmp[3:6]]
box[2,:] = [float(x) for x in tmp[6:9]]

print (box)

AtomList = []
xyz = np.zeros(shape=(natom,3))
for k in range(natom):
    tmp = f.readline()
    tmp = tmp.split()
    AtomList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
f.close


def write_POSCAR(fnamePOSCAR, box, xyz, AtomList):
    f2 = open(fnamePOSCAR, "w")
    f2.write("%1s\n" %( "COMMENT" ))
    f2.write("%15.9f\n" %( 1.0 ))
    f2.write("%15.9f %15.9f %15.9f\n" %( box[0,0], box[0,1], box[0,2] ))
    f2.write("%15.9f %15.9f %15.9f\n" %( box[1,0], box[1,1], box[1,2] ))
    f2.write("%15.9f %15.9f %15.9f\n" %( box[2,0], box[2,1], box[2,2] ))
    
    syms, counts_syms = np.unique(AtomList, return_counts=True)
    
    print (syms)
    print (counts_syms)
    
    f3 = open("ntype.dat", "w")
    for item in syms:
        f2.write("%s %s" %(" ", item))
        f3.write("%s\n" %( item))
    
    f2.write("\n")
    for item in counts_syms:
        f2.write("%s %s" %(" ", item))
    f2.write("\n")
    f2.write("%1s\n" %( "Direct" ))
    
    pxyz = np.dot(xyz, np.linalg.inv(box))
    
    for item in syms:
        for k in range(natom):
            if (item == AtomList[k]):
                f2.write("%15.9f %15.9f %15.9f\n" %( pxyz[k,0], pxyz[k,1], pxyz[k,2] ))
    f2.close()

def generate_xyz_POSCAR(fnamePOSCAR, OH_bond, box, xyz, AtomList):
    nxyz = xyz
    nmole = natom/3
    for imole in range(nmole):
        idO = 3*imole + 2
        idH1= 3*imole
        idH2= 3*imole + 1
        xyz_O  = nxyz[idO ,:]
        xyz_H1 = nxyz[idH1,:]
        xyz_H2 = nxyz[idH2,:]
        v1 = xyz_H1 - xyz_O
        d1 = np.sqrt(np.dot(v1,v1))
        v2 = xyz_H2 - xyz_O
        d2 = np.sqrt(np.dot(v2,v2))
        nxyz[idH1,:] = xyz_O + v1/d1*OH_bond
        nxyz[idH2,:] = xyz_O + v2/d2*OH_bond
        print (imole, idO, idH1, idH2)
    write_POSCAR(fnamePOSCAR, box, nxyz, AtomList)

fnamePOSCAR = "POSCAR"
generate_xyz_POSCAR(fnamePOSCAR, OH_bond, box, xyz, AtomList)

