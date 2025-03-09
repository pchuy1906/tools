import os
import numpy as np
from numpy import *
from numpy import matrix
from numpy import linalg
import commands

import argparse
parser = argparse.ArgumentParser(description='export FORCES, ENERGIES, and STRESS of file xyzfes')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')

args        = parser.parse_args()
file_xyz     = args.file_xyz

print ("")
print ("Orthohombic cell case")
print ("")
print ("Working only for constant number of atom")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"r")

status, nline = commands.getstatusoutput(" wc " + file_xyz + " | awk '{print $1}' ")
nline = int(nline)

nstruc = 10000
istruc = 0

Ha2kcalmol=627.509
Bohr2Angstrom=0.529177

f2 = open("reference_force_stress_energy.dat", "w")

all_force = []

while (istruc < nstruc):
    #print ("istruc, nstruc = ", istruc, nstruc)
    natom = f.readline()
    natom = int(natom)
    
    nstruc = nline / (natom + 2)
    
    #print ("the number of atom:", natom)
    
    box = np.zeros(shape=(3,3))
    force = np.zeros(shape=(3*natom))
    stress = np.zeros(shape=(6))
    energy = np.zeros(shape=(3))

    tmp = f.readline()
    tmp_ = tmp.split()

    for k in range(3):
        box[k,k] = float(tmp_[k])

    if (len(tmp_)==4):
	for k in range(3):
	    energy[k] = float(tmp_[3])
    elif (len(tmp_)==10):
        for k in range(3):
            energy[k] = float(tmp_[9])
        for k in range(6):
	    stress[k] = float(tmp_[3+k])

    myList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline()
        tmp = tmp.split()
        myList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
        force[3*k+0] = float(tmp[4])* Ha2kcalmol/Bohr2Angstrom
        force[3*k+1] = float(tmp[5])* Ha2kcalmol/Bohr2Angstrom
        force[3*k+2] = float(tmp[6])* Ha2kcalmol/Bohr2Angstrom
  	# VASP has force (eV/A), while ChIMES used (Ha/Bohr) for input; and (kcal/mol/A) for output
        f2.write("%15.9f %1s\n" %( force[3*k+0], " FORCE" ))
        f2.write("%15.9f %1s\n" %( force[3*k+1], " FORCE" ))
        f2.write("%15.9f %1s\n" %( force[3*k+2], " FORCE" ))
    # export the stress tensors
    if (len(tmp_)==10):
        f2.write("%15.9f %1s\n" %( stress[0], " STRESS" ))
        f2.write("%15.9f %1s\n" %( stress[3], " STRESS" ))
        f2.write("%15.9f %1s\n" %( stress[5], " STRESS" ))

        f2.write("%15.9f %1s\n" %( stress[3], " STRESS" ))
        f2.write("%15.9f %1s\n" %( stress[1], " STRESS" ))
        f2.write("%15.9f %1s\n" %( stress[4], " STRESS" ))

        f2.write("%15.9f %1s\n" %( stress[5], " STRESS" ))
        f2.write("%15.9f %1s\n" %( stress[4], " STRESS" ))
        f2.write("%15.9f %1s\n" %( stress[2], " STRESS" ))

        f2.write("%15.9f %1s\n" %( energy[0], " ENERGY" ))
        f2.write("%15.9f %1s\n" %( energy[1], " ENERGY" ))
        f2.write("%15.9f %1s\n" %( energy[2], " ENERGY" ))

    if (len(tmp_)==4):
        f2.write("%15.9f %1s\n" %( energy[0], " ENERGY" ))
        f2.write("%15.9f %1s\n" %( energy[1], " ENERGY" ))
        f2.write("%15.9f %1s\n" %( energy[2], " ENERGY" ))

    all_force = np.concatenate([all_force,force])
    istruc += 1
f.close
print 
print ("the number of structures:", nstruc)
print ("the number of force components:", len(all_force))
print ("min-force", min(all_force), " kcal/mol-A")
print ("max-force", max(all_force), " kcal/mol-A")
print 

bins_inp = np.linspace( -4000, 4000, num=200)
#from matplotlib import pyplot as plt 
#plt.hist(all_force, bins = bins) 
#plt.show()
hist,bins = np.histogram(all_force, bins = bins_inp) 
res = np.vstack((bins[:-1], hist)).T

np.savetxt('hist.dat', res, fmt='%15.9f %15.9f')

