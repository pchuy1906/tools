import os
import numpy as np
import math

import argparse

parser = argparse.ArgumentParser(description='xyz_2_gen')

# Arguments supported by the code.
parser.add_argument("--file_xyz",      default='file.xyz', help='file.xyz')
parser.add_argument("--option_cell",   default='cell_3',   help='NON_ORTHO/cell_3/cell_9')
args        = parser.parse_args()
file_xyz     = args.file_xyz
option_cell  = args.option_cell

f  = open(file_xyz,      "r")
f2 = open("xyz2gen.gen", "w")

nstep = 0
while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    tmp = tmp.split()
    nstep += 1
    natom = int(tmp[0])
    #print ("number of atom ", natom)
    
    tmp = f.readline().split()
    unitcell = np.zeros(shape=(9))
    if (option_cell=="cell_9"):
        unitcell = tmp[0:9]
    elif (option_cell=="cell_3"):
        unitcell[0] = float(tmp[0])
        unitcell[4] = float(tmp[1])
        unitcell[8] = float(tmp[2])
    else:
        print ("unknown option")

    unitcell = np.array(unitcell)
    unitcell = unitcell.reshape((3, 3))

    xyz = np.zeros(shape=(natom,3))
    atomName = []
    for k in range(natom):
        tmp = f.readline().split()
        atomName.append(tmp[0])
        xyz[k,:] =  tmp[1:4]
    
    syms, counts_syms = np.unique(atomName, return_counts=True)
    
    #print (syms)
    asym_list = ' '.join(syms)
    #print (asym_list)
    #print (counts_syms)
    
    f2.write("%-d %4s\n" %(natom, "S"))
    f2.write("%-s \n" %(asym_list))
    
    for k in range(natom):
        ind = np.where(syms == atomName[k])[0][0]+1
        f2.write("%-5d %5d %15.9f %15.9f %15.9f\n" %(k+1, ind, xyz[k,0], xyz[k,1], xyz[k,2] ))
    
    f2.write("%15.9f%15.9f%15.9f\n" % (0.0,0.0,0.0))
    for ixyz in range(3):
        for icellxyz in unitcell[ixyz,:]:
            f2.write("%15.9f" % icellxyz)
        f2.write("\n")
    
f.close()
print ("number of configurations", nstep)
