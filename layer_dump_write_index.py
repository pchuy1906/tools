import numpy as np

'''
structure of the dump is as follow:
-----
nH_host
-----
2 n_wat
-----
nO_host
-----
n_wat
-----
nSi_host
'''
 
# read input 
import argparse
parser = argparse.ArgumentParser(description='split file xyz')
# Arguments supported by the code.
parser.add_argument("--file_dump",       default='file.dump', help='file with format dump lammps')
parser.add_argument("--fPOSCAR",         default='POSCAR',    help='file POSCAR for whole system')
parser.add_argument("--nhost", type=int, default=1,           help='number of host atoms')
parser.add_argument("--nwat",  type=int, default=1,           help='number of water molecules')
parser.add_argument('--zboundary', nargs='+', type=float)


args    = parser.parse_args()
file_dump = args.file_dump
fPOSCAR   = args.fPOSCAR
nhost     = args.nhost
nwat      = args.nwat
zboundary = args.zboundary

f  = open(file_dump ,"rt")
nconf = 0

f2 = open(fPOSCAR)
lines = f2.readlines()
# line number 7 is: nH, nO, nSi
nH_all, nO_all, nSi_all = lines[6].split()
nH_all = int(nH_all)
nO_all = int(nO_all)
nSi_all = int(nSi_all)
print (nH_all, nO_all, nSi_all)

nH_wat = 2*nwat
nO_wat = 1*nwat
nH_host = nH_all - nH_wat
nO_host = nO_all - nO_wat
nSi_host = nSi_all

natom_real = 3*nwat

nlayer = len(zboundary)-1
natom_layer  = np.zeros(shape=(nlayer), dtype=int)

while (nconf<1):
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    tmp  = f.readline()
    tmp  = f.readline()
    tmp  = f.readline()
    natom = int(tmp)

    tmp = f.readline()
    tmp = f.readline()
    tmp = f.readline()
    tmp = f.readline()

    tmp = f.readline().split()
    tmp = np.array(tmp)
    idx  = np.where(tmp=='x')[0][0]
    idy  = np.where(tmp=='y')[0][0]
    idz  = np.where(tmp=='z')[0][0]
    idi  = np.where(tmp=='type')[0][0]
    idvx  = np.where(tmp=='vx')[0][0]
    idvy  = np.where(tmp=='vy')[0][0]
    idvz  = np.where(tmp=='vz')[0][0]

    id_atoms  = np.zeros(shape=(natom), dtype=int)

    x  = np.zeros(shape=(natom))
    y  = np.zeros(shape=(natom))
    z  = np.zeros(shape=(natom))
    vx  = np.zeros(shape=(natom))
    vy  = np.zeros(shape=(natom))
    vz  = np.zeros(shape=(natom))
    iid  = np.zeros(shape=(natom), dtype=int)
    for i in range(natom):
        tmp = f.readline().split()
        idt = int(tmp[0])
        x[idt-1]  = float(tmp[idx-2])
        y[idt-1]  = float(tmp[idy-2])
        z[idt-1]  = float(tmp[idz-2])
        vx[idt-1]  = float(tmp[idvx-2])
        vy[idt-1]  = float(tmp[idvy-2])
        vz[idt-1]  = float(tmp[idvz-2])
        iid[idt-1] =  int(tmp[idi-2])

    ncount = 0

    
    # nH_host  atoms: do nothing
    # nH_water atoms: wait for oxygen
    # nO_host  atoms: do nothing
    # nO_water atoms:      
    index_begin = nH_host + 2*nwat + nO_host
    index_end   = nH_host + 2*nwat + nO_host + nwat
    for i in range(index_begin, index_end):
        ncount += 1
        for j in range(len(zboundary)-1):
            if zboundary[j] < z[i] and z[i] <= zboundary[j+1]:
                id_atoms[i] = j+1
                natom_layer[j] +=3
    nconf += 1
f.close()

np.savetxt('atom_index', id_atoms, fmt='%d')
np.savetxt('natom_index',natom_layer, fmt='%d')


print ("number of frame ", nconf)
