import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import argparse
parser = argparse.ArgumentParser(description='modify XYZ for DFTB set up')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')

args        = parser.parse_args()
file_xyz           = args.file_xyz

print ("-----------------------------------------------------------")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"rt")

istruc = 0
while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break
    #v, e = [float(x) for x in line.split()[:2]]

    natom = int(tmp)
    
    box = np.zeros(shape=(3,3))
    force = np.zeros(shape=(3*natom))
    stress = np.zeros(shape=(6))
    energy = np.zeros(shape=(3))

    tmp = f.readline().split()
    print (tmp)
    
    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline().split()
        atomList.append(tmp[0])
        xyz[k,0],xyz[k,1],xyz[k,2] =  float(tmp[1]),float(tmp[2]),float(tmp[3])
        force[3*k+0],force[3*k+1],force[3*k+2] = float(tmp[4]),float(tmp[5]),float(tmp[6])
    istruc += 1
f.close

print ("The number of structure is %d" %istruc)


#if (not fit_AE_only):
#
#    f2 = open("substract_atomic_energy_"+file_xyz, "w")
#    
#    f  = open(file_xyz ,"rt")
#    istruc=0
#    
#    while True:
#        tmp  = f.readline()
#        line = tmp.strip()
#        if line == '': break
#    
#        natom = int(tmp)
#        f2.write(tmp)
#    
#        tmp = f.readline()
#        tmp_ = tmp.split()
#        energy = (float(tmp_[3])-np.dot(Amatrix[istruc,:],Atomic_energies[:]))
#        tmp_= np.array(tmp_)
#        f2.write("%12.3f %12.3f %12.3f %15.9f\n" %(float(tmp_[0]), float(tmp_[1]), float(tmp_[2]), energy))
#    
#        for k in range(0,natom):
#            tmp = f.readline()
#            f2.write(tmp)
#        istruc+=1
#    f.close
    
