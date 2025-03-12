import numpy as np

import argparse
parser = argparse.ArgumentParser(description='export atomic forces, energies, and stresses of file xyzfes')

# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file_XYZ format xyz')

args        = parser.parse_args()
file_xyz     = args.file_xyz

print ("")
print ("Export atomic forces, energies, and stresses")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"rt")

Ha2kcalmol=627.509
Bohr2Angstrom=0.529177

f2 = open("reference_force_stress_energy.dat", "w")
f3 = open("number_of_atom.dat", "w")

all_force = []

istruc = 0
while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp)
    
    box = np.zeros(shape=(3,3))
    force = np.zeros(shape=(3*natom))
    stress = np.zeros(shape=(6))
    energy = np.zeros(shape=(3))

    tmp = f.readline().split()

    ncomment = len(tmp)
    # the last column is the energy
    energy = float(tmp[ncomment-1])

    istruc += 1

    myList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline().split()
        myList.append(tmp[0])
        xyz[k,0], xyz[k,1], xyz[k,2] =  float(tmp[1]), float(tmp[2]), float(tmp[3])
        force[3*k+0] = float(tmp[4])* Ha2kcalmol/Bohr2Angstrom
        force[3*k+1] = float(tmp[5])* Ha2kcalmol/Bohr2Angstrom
        force[3*k+2] = float(tmp[6])* Ha2kcalmol/Bohr2Angstrom
  	# VASP has force (eV/A), while ChIMES used (Ha/Bohr) for input; and (kcal/mol/A) for output
        TMP_STR = str(istruc)+"_FORCE_"+tmp[0]+"_x"
        f2.write("%15.9f %1s\n" %( force[3*k+0], TMP_STR ))
        TMP_STR = str(istruc)+"_FORCE_"+tmp[0]+"_y"
        f2.write("%15.9f %1s\n" %( force[3*k+1], TMP_STR ))
        TMP_STR = str(istruc)+"_FORCE_"+tmp[0]+"_z"
        f2.write("%15.9f %1s\n" %( force[3*k+2], TMP_STR ))
        #f3.write("%15.9f\n" %( 1.0 ))
        #f3.write("%15.9f\n" %( 1.0 ))
        #f3.write("%15.9f\n" %( 1.0 ))

    # export the energy
    TMP_STR = str(istruc)+"_ENERGY"
    f2.write("%15.9f %1s\n" %( energy, TMP_STR ))
    f3.write("%12.0f\n" %( natom ))

    all_force = np.concatenate([all_force,force])
f.close
print 
print ("the number of structures:", istruc)
print ("the number of force components:", len(all_force))
print ("min-force", min(all_force), " kcal/mol-A")
print ("max-force", max(all_force), " kcal/mol-A")
print 

#bins_inp = np.linspace( -4000, 4000, num=200)
#from matplotlib import pyplot as plt 
#plt.hist(all_force, bins = bins) 
#plt.show()
#hist,bins = np.histogram(all_force, bins = bins_inp) 
#res = np.vstack((bins[:-1], hist)).T
#
#np.savetxt('hist.dat', res, fmt='%15.9f %15.9f')
