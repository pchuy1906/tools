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
parser.add_argument("--weight_force",  type=float, default=1.0, help='weight forces')
parser.add_argument("--weight_energy", type=float, default=300.0, help='weight energy')
parser.add_argument("--weight_stress", type=float, default=500.0, help='weight stress')


args        = parser.parse_args()
file_xyz     = args.file_xyz
weight_force     = args.weight_force
weight_energy    = args.weight_energy
weight_stress    = args.weight_stress

print ("")
print ("The general case")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"rt")

Ha2kcalmol=627.509
Bohr2Angstrom=0.529177

f2 = open("new_weight.dat", "w")

all_force = []

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

    ioption = 3

    if (ioption==1):
        box[0,:] = [float(x) for x in tmp[1:4]]
        box[1,:] = [float(x) for x in tmp[4:7]]
        box[2,:] = [float(x) for x in tmp[7:10]]

        stress = [float(x) for x in tmp[10:16]]
        #print (stress)
        for k in range(3):
            energy[k] = float(tmp[16])
        #print (energy)
        fitStress = True
    elif (ioption == 2):
        box[0,:] = [float(x) for x in tmp[1:4]]
        box[1,:] = [float(x) for x in tmp[4:7]]
        box[2,:] = [float(x) for x in tmp[7:10]]

        for k in range(3):
            energy[k] = float(tmp[10])
        #print (energy)
        fitStress = False
    elif (ioption==3):
        box[0,0] = float(tmp[0])
        box[1,1] = float(tmp[1])
        box[2,2] = float(tmp[2])

        try:
            stress = [float(x) for x in tmp[3:9]]
            #print (stress)
            for k in range(3):
                energy[k] = float(tmp[9])
            #print (energy)
            fitStress = True
        except:
            fitStress = False
            for k in range(3):
                energy[k] = float(tmp[3])


    #print (box)


    myList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline().split()
        myList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
        force[3*k+0] = float(tmp[4])* Ha2kcalmol/Bohr2Angstrom
        force[3*k+1] = float(tmp[5])* Ha2kcalmol/Bohr2Angstrom
        force[3*k+2] = float(tmp[6])* Ha2kcalmol/Bohr2Angstrom
  	# VASP has force (eV/A), while ChIMES used (Ha/Bohr) for input; and (kcal/mol/A) for output
        f2.write("%15.9f\n" %( weight_force ))
        f2.write("%15.9f\n" %( weight_force ))
        f2.write("%15.9f\n" %( weight_force ))
    # export the stress tensors
    if fitStress:
        f2.write("%15.9f \n" %( weight_stress ))
        f2.write("%15.9f \n" %( weight_stress ))
        f2.write("%15.9f \n" %( weight_stress ))
    
        f2.write("%15.9f \n" %( weight_stress ))
        f2.write("%15.9f \n" %( weight_stress ))
        f2.write("%15.9f \n" %( weight_stress ))
    
        f2.write("%15.9f \n" %( weight_stress ))
        f2.write("%15.9f \n" %( weight_stress ))
        f2.write("%15.9f \n" %( weight_stress ))
    
    # export energy
    f2.write("%15.9f \n" %( weight_energy/float(natom) ))
    f2.write("%15.9f \n" %( weight_energy/float(natom) ))
    f2.write("%15.9f \n" %( weight_energy/float(natom) ))

    all_force = np.concatenate([all_force,force])
    istruc += 1
f.close
print 
print ("the number of structures:", istruc)
print ("the number of force components:", len(all_force))
print ("min-force", min(all_force), " kcal/mol-A")
print ("max-force", max(all_force), " kcal/mol-A")
print 


f3 = open("force_min_max.dat", "w")
if (min(all_force) < -5000.0):
    f3.write("%s \n" %("NO"))
if (max(all_force) >  5000.0):
    f3.write("%s \n" %("NO"))
f3.close

#bins_inp = np.linspace( -4000, 4000, num=200)
#from matplotlib import pyplot as plt 
#plt.hist(all_force, bins = bins) 
#plt.show()
#hist,bins = np.histogram(all_force, bins = bins_inp) 
#res = np.vstack((bins[:-1], hist)).T
#
#np.savetxt('hist.dat', res, fmt='%15.9f %15.9f')

