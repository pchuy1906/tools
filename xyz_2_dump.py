import os
import numpy as np
import math

import argparse

parser = argparse.ArgumentParser(description='xyz_2_dump')

# Arguments supported by the code.
parser.add_argument("--file_xyz",      default='file.xyz', help='file.xyz')
parser.add_argument("--option_cell",   default='cell_3',   help='NON_ORTHO/cell_3/cell_9')
args        = parser.parse_args()
file_xyz     = args.file_xyz
option_cell  = args.option_cell

f  = open(file_xyz,      "r")
f2 = open("lammps.dump", "w")

nstep = 0
while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    tmp = tmp.split()
    nstep += 100
    natom = int(tmp[0])
    #print ("number of atom ", natom)
    
    tmp = f.readline().split()
    unitcell = np.zeros(shape=(9))
    if (option_cell=="cell_9"):
        unitcell = tmp[0:9]
        unitcell = np.array(unitcell)
        unitcell = unitcell.reshape((3, 3))

    xyz = np.zeros(shape=(natom,3))
    atomName = []
    for k in range(natom):
        tmp = f.readline().split()
        atomName.append(tmp[0])
        xyz[k,:] =  tmp[1:4]
    
    f2.write("ITEM: TIMESTEP\n")
    f2.write("%-d\n" %(nstep))

    f2.write("ITEM: NUMBER OF ATOMS\n")
    f2.write("%-d\n" %(natom))

    f2.write("ITEM: BOX BOUNDS pp pp pp\n")
    for ixyz in range(3):
        f2.write("%15.9f %s\n" % (0.0, unitcell[ixyz,ixyz]) )

    f2.write("ITEM: ATOMS id type xu yu zu\n")
    for i in range(natom):
        f2.write("%8d" % (i+1) )
        f2.write("%4s" % atomName[i] )
        for ixyz in range(3):
            f2.write("%15.9f" % xyz[i,ixyz] )
        f2.write("\n")


f.close()
