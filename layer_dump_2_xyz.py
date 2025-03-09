import numpy as np
import time
import datetime
now = datetime.datetime.now()

print (" starting time:", now)

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
parser.add_argument("--file_dump",        default='file.dump',   help='file with format dump lammps')
parser.add_argument("--fPOSCAR",          default='POSCAR',      help='file POSCAR for whole system')
parser.add_argument("--nhost", type=int,  default=1,             help='number of host atoms')
parser.add_argument("--nwat",  type=int,  default=1,             help='number of water molecules')
parser.add_argument("--file_atom_index",  default='atom_index',  help='file with atom indexes')
parser.add_argument("--file_natom_index", default='natom_index', help='file with natom for each layer')

args    = parser.parse_args()
file_dump        = args.file_dump
fPOSCAR          = args.fPOSCAR
nhost            = args.nhost
nwat             = args.nwat
file_atom_index  = args.file_atom_index
file_natom_index = args.file_natom_index

f  = open(file_dump ,"rt")
nconf = 0

atom_index = np.loadtxt(file_atom_index, dtype=int)
natom_index = np.loadtxt(file_natom_index, dtype=int)

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

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    tmp  = f.readline()
    tmp  = f.readline()
    tmp  = f.readline()
    natom = int(tmp)

    for i in range(len(natom_index)):
        fname = "layer_" + str(i+1) + ".xyz"
        f2 = open(fname, "a")
        f2.write("%6d\n" %( natom_index[i] ))
        f2.close()

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


    # print Oxygen first, with type 1
    index_begin = nH_host + 2*nwat + nO_host
    index_end   = nH_host + 2*nwat + nO_host + nwat
    for i in range(index_begin, index_end):
        ncount += 1

        fname = "layer_" + str(atom_index[i]) + ".xyz"
        f2 = open(fname, "a")
        f2.write("%s %15.6f %15.6f %15.6f\n" %( "O ", x[i], y[i], z[i] ))
        j = nH_host + 2*(i - index_begin) + 0
        f2.write("%s %15.6f %15.6f %15.6f\n" %( "H ", x[j], y[j], z[j] ))
        j = nH_host + 2*(i - index_begin) + 1
        f2.write("%s %15.6f %15.6f %15.6f\n" %( "H ", x[j], y[j], z[j] ))
        f2.close()

    nconf += 1
f.close()

print ("number of frame ", nconf)

now = datetime.datetime.now()
print (" end time:", now)


