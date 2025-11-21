import numpy as np
 
# read input 
import argparse
parser = argparse.ArgumentParser(description='split file xyz')
# Arguments supported by the code.
parser.add_argument("--file_dump",   default='file.dump', help='file with format dump lammps')
parser.add_argument("--read_charge", default=False, action="store_true", help="read point charge")

args    = parser.parse_args()
file_dump     = args.file_dump
read_charge   = args.read_charge

atm_2_GPa = 0.000101325
bar_2_GPa = 0.0001

A_2_Bohr = 1.889725989
kcalmol_2_Ha = 0.00159362
kcalmol_A_2_Ha_Bohr = kcalmol_2_Ha / A_2_Bohr

eV_2_Ha = 0.0367502 
eV_A_2_Ha_Bohr = eV_2_Ha / A_2_Bohr

eV_2_kcalmol = 23.0609 


f  = open(file_dump ,"rt")
nconf = 1

fname = "charges.dat"
f2 = open( fname, "w")

while True:

    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break
    tmp  = f.readline()
    tmp  = f.readline()
    tmp  = f.readline()
    natom = int(tmp)

    tmp = f.readline()
    if "xy" in tmp:
        #print ("NON_ORTHO")

        tmp = f.readline().split()
        [xlo_bound, xhi_bound, xy] = [float(x) for x in tmp]

        tmp = f.readline().split()
        [ylo_bound, yhi_bound, xz] = [float(x) for x in tmp]

        tmp = f.readline().split()
        [zlo_bound, zhi_bound, yz] = [float(x) for x in tmp]

        zlo = zlo_bound
        zhi = zhi_bound
        ylo = ylo_bound - min(0.0,yz)
        yhi = yhi_bound - max(0.0,yz)
        xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
        xhi = xhi_bound - max(0.0,xy,xz,xy+xz)

        lx = xhi-xlo
        ly = yhi-ylo
        lz = zhi-zlo
    else:
        tmp = f.readline().split()
        #print (tmp)
        lx = float(tmp[1])-float(tmp[0])
        tmp = f.readline().split()
        #print (tmp)
        ly = float(tmp[1])-float(tmp[0])
        tmp = f.readline().split()
        #print (tmp)
        lz = float(tmp[1])-float(tmp[0])
        #print (lx,ly,lz)
        xy = 0.0
        xz = 0.0
        yz = 0.0

    cell_9 = [lx,0,0, xy,ly,0, xz,yz,lz]

    tmp = f.readline().split()
    tmp = np.array(tmp)
    if (read_charge):
        idq  = np.where(tmp=='q')[0][0]
    #idfx = np.where(tmp=='fx')[0][0]
    #idfy = np.where(tmp=='fy')[0][0]
    #idfz = np.where(tmp=='fz')[0][0]
    idx  = np.where(tmp=='x')[0][0]
    idy  = np.where(tmp=='y')[0][0]
    idz  = np.where(tmp=='z')[0][0]
    idi  = np.where(tmp=='type')[0][0]

    #print (idq, idfx, idfy, idfz, idx, idy, idz, idi)
    q  = np.zeros(shape=(natom))
    x  = np.zeros(shape=(natom))
    y  = np.zeros(shape=(natom))
    z  = np.zeros(shape=(natom))
    #fx = np.zeros(shape=(natom))
    #fy = np.zeros(shape=(natom))
    #fz = np.zeros(shape=(natom))
    strs = ["XXXX" for x in range(natom)]
    atype = np.array(strs)
    for i in range(natom):
        tmp = f.readline().split()
        idt = int(tmp[0])
        if (read_charge):
            q[idt-1]  = float(tmp[idq-2])
        x[idt-1]  = float(tmp[idx-2])
        y[idt-1]  = float(tmp[idy-2])
        z[idt-1]  = float(tmp[idz-2])

        f_convert = 1.0
        p_convert = 1.0
        e_convert = 1.0

        #fx[idt-1] = float(tmp[idfx-2]) * f_convert
        #fy[idt-1] = float(tmp[idfy-2]) * f_convert
        #fz[idt-1] = float(tmp[idfz-2]) * f_convert
        iid = int(tmp[idi-2])

    if (read_charge):
        f2.write("%d " %(nconf))
        for i in range(natom):
            f2.write("%6.3f " %(q[i]))
        f2.write("\n")

    nconf += 1


f2.close()


