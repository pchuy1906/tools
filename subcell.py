import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg

import argparse
parser = argparse.ArgumentParser(description='making the larger system')
# Arguments supported by the code.
parser.add_argument("--file_xyz",     default='file.xyz', help='file_XYZ format xyz')
parser.add_argument("--cell_type",    default='cell_3',   help='cell_3/cell_9/NON_ORTHO')
parser.add_argument("--Nx", type=int, default=2,          help='Nx')
parser.add_argument("--Ny", type=int, default=2,          help='Ny')
parser.add_argument("--Nz", type=int, default=2,          help='Nz')

args        = parser.parse_args()
file_xyz     = args.file_xyz
cell_type    = args.cell_type
Nx           = args.Nx
Ny           = args.Ny
Nz           = args.Nz

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

if (cell_type=="cell_3"):
    print ("")
    print ("Orthohombic cell case")
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
elif (cell_type=="cell_9"):
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
    box = [float(x) for x in tmp[:9]]
    box = np.array(box)
    box = box.reshape((3, 3))

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
elif (cell_type=="NON_ORTHO"):
    print ("")
    print ("NON_ORTHO")
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
    box = [float(x) for x in tmp[1:10]]
    box = np.array(box)
    box = box.reshape((3, 3))

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
else:
    print ("Not implemented yet!")
    exit()

newbox = np.zeros(shape=(3,3))
newbox[0,:] = box[0,:]*Nx
newbox[1,:] = box[1,:]*Ny
newbox[2,:] = box[2,:]*Nz

np.savetxt('newbox', newbox, fmt='%15.9f %15.9f %15.9f')

f2 = open("bigger_"+str(Nx)+"_"+str(Ny)+"_"+str(Nz)+".xyz", "w")
f2.write("%4d\n" %( natom*Nx*Ny*Nz ))
if (cell_type=='cell_3'):
    f2.write("%15.9f %15.9f %15.9f\n" %( newbox[0,0], newbox[1,1], newbox[2,2] ))
else:
    f2.write("%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n" %( newbox[0,0], newbox[0,1], newbox[0,2], newbox[1,0], newbox[1,1], newbox[1,2], newbox[2,0], newbox[2,1], newbox[2,2] ))

for k1 in range(1,Nx+1):
    for k2 in range(1,Ny+1):
        for k3 in range(1,Nz+1):
            trans = float(k1-1)*box[0,:] + float(k2-1)*box[1,:] + float(k3-1)*box[2,:]
            for k4 in range(0,natom):
                tmpp  = xyz[k4,:]+trans[:]
                f2.write("%4s %15.9f %15.9f %15.9f\n" %( myList[k4], tmpp[0], tmpp[1], tmpp[2] ))
f2.close()
