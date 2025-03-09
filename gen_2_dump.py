import os
import numpy as np
import math

import argparse

parser = argparse.ArgumentParser(description='gen_2_dump')

# Arguments supported by the code.
parser.add_argument("--file_input", default='file.gen', help='gen_2_xyz')
parser.add_argument("--cell_out",   default='cell_3',   help='NON_ORTHO/cell_3/cell_9')
args        = parser.parse_args()
file_input   = args.file_input
cell_out     = args.cell_out

f  = open(file_input ,   "r")
f2 = open("lammps.dump", "w")

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    tmp = tmp.split()
    nstep = int(tmp[3])
    natom = int(tmp[0])
    #print ("number of atom ", natom)
    
    atomNameList = f.readline().split()
    #print (atomNameList)
    
    xyz = np.zeros(shape=(natom,3))
    atomName = []
    for k in range(natom):
        tmp = f.readline().split()
        ID  = int(tmp[1])-1
        atomName.append(atomNameList[ID])
        xyz[k,:] =  tmp[2:5]
    
    #print (xyz)
    
    tmp = f.readline()
    
    unitcell = np.zeros(shape=(3,3))
    tmp = f.readline().split()
    unitcell[0,:] = tmp
    tmp = f.readline().split()
    unitcell[1,:] = tmp
    tmp = f.readline().split()
    unitcell[2,:] = tmp
    f2.write("ITEM: TIMESTEP\n")
    f2.write("%-d\n" %(nstep))

    f2.write("ITEM: NUMBER OF ATOMS\n")
    f2.write("%-d\n" %(natom))

    f2.write("ITEM: BOX BOUNDS pp pp pp\n")
    for ixyz in range(3):
        f2.write("%15.9f %15.9f\n" % (0.0, unitcell[ixyz,ixyz]) )

    f2.write("ITEM: ATOMS id type xu yu zu\n")
    for i in range(natom):
        f2.write("%8d" % (i+1) )
        f2.write("%4s" % atomName[i] )
        for ixyz in range(3):
            f2.write("%15.9f" % xyz[i,ixyz] )
        f2.write("\n")


f.close()
