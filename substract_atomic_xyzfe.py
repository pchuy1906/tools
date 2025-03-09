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
parser.add_argument("--file_atomic_energy", default='atomic_energies.dat', help='input/output the atomic energies')
parser.add_argument("--fit_atomic_energy", action='store_true', help='fit or use atomic energies')
parser.add_argument("--fit_AE_only", action='store_true', help='fit atomic energies ONLY')

args        = parser.parse_args()
file_xyz           = args.file_xyz
file_atomic_energy = args.file_atomic_energy
fit_atomic_energy  = args.fit_atomic_energy
fit_AE_only        = args.fit_AE_only

print ("-----------------------------------------------------------")
print ("")
print ("Orthohombic cell case only")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"rt")

def number_CHNO(atom_order, molecule_order, molecule_natom):
    res = np.array([0,0,0,0])
    for id_tmp in range(len(molecule_order)):
        atom_loc = np.where(atom_order == molecule_order[id_tmp])
        res[atom_loc] = molecule_natom[id_tmp]
    return res


istruc = 0
atom_order = np.array(['C','H','N','O'])

#print (number_CHNO(atom_order=atom_order, molecule_order=values, molecule_natom=counts))
Amatrix = np.array([])
bmatrix = []

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

    tmp = f.readline()
    tmp_ = tmp.split()

    for k in range(3):
        box[k,k] = float(tmp_[k])
	energy[k] = float(tmp_[3])/float(natom)
    bmatrix.append(float(tmp_[3]))

    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline()
        tmp = tmp.split()
        atomList.append(tmp[0])
        xyz[k,0],xyz[k,1],xyz[k,2] =  float(tmp[1]),float(tmp[2]),float(tmp[3])
        force[3*k+0],force[3*k+1],force[3*k+2] = float(tmp[4]),float(tmp[5]),float(tmp[6])
    molecule_order, molecule_natom = np.unique(atomList, return_counts=True)
    Amatrix = np.append(Amatrix, number_CHNO(atom_order, molecule_order, molecule_natom) )
    istruc += 1
f.close

Amatrix = Amatrix.reshape((istruc, 4))
bmatrix = np.array(bmatrix)

print 
print ("the number of structures:", istruc)
print 


if (fit_atomic_energy):
    Atomic_energies = np.linalg.lstsq(Amatrix, bmatrix,rcond=None)[0]
    print ("\nThe calculated atomics energies are:")
    print (Atomic_energies)
    np.savetxt(file_atomic_energy, Atomic_energies)
else:
    print ("\nRead the atomic energies:")
    Atomic_energies = np.loadtxt(file_atomic_energy)
    print (Atomic_energies)


if (not fit_AE_only):

    f2 = open("substract_atomic_energy_"+file_xyz, "w")
    
    f  = open(file_xyz ,"rt")
    istruc=0
    
    while True:
        tmp  = f.readline()
        line = tmp.strip()
        if line == '': break
    
        natom = int(tmp)
        f2.write(tmp)
    
        tmp = f.readline()
        tmp_ = tmp.split()
        energy = (float(tmp_[3])-np.dot(Amatrix[istruc,:],Atomic_energies[:]))
        tmp_= np.array(tmp_)
        f2.write("%12.3f %12.3f %12.3f %15.9f\n" %(float(tmp_[0]), float(tmp_[1]), float(tmp_[2]), energy))
    
        for k in range(0,natom):
            tmp = f.readline()
            f2.write(tmp)
        istruc+=1
    f.close
    
