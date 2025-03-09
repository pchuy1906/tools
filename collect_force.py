import numpy as np
 
# read input 
import argparse
parser = argparse.ArgumentParser(description='split file xyz')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file with format xyzf, xyzfe, xyzfes')
args    = parser.parse_args()
file_xyz          = args.file_xyz

f  = open(file_xyz ,"rt")
nconf = 1

Ha2kcalmol=627.509
Bohr2Angstrom=0.529177


all_fxyz = np.array([])

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp)

    tmp = f.readline()

    fx = np.zeros(shape=(natom))
    fy = np.zeros(shape=(natom))
    fz = np.zeros(shape=(natom))

    for i in range(natom):
        tmp = f.readline().split()
        fx[i] = float(tmp[4])
        fy[i] = float(tmp[5])
        fz[i] = float(tmp[6])

    fxyz = np.concatenate((fx, fy, fz))
 
    tmp_fxyz = np.concatenate((all_fxyz,fxyz))
    all_fxyz = tmp_fxyz

    nconf += 1

print (len(all_fxyz))
all_fxyz = all_fxyz * Ha2kcalmol/Bohr2Angstrom

np.savetxt('DFT_forces.dat', all_fxyz, fmt='%15.9f')

