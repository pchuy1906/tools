import numpy as np
 
# read input 
import argparse
parser = argparse.ArgumentParser(description='split file xyz')
# Arguments supported by the code.
parser.add_argument("--file_dump",   default='file.dump', help='file with format dump lammps')
parser.add_argument("--cell_type",   default='NON_ORTHO', help='cell_type (cell_3/cell_9/NON_ORTHO)')
parser.add_argument('--atom_types', nargs='+')

args    = parser.parse_args()
file_dump     = args.file_dump
atom_types    = args.atom_types
cell_type     = args.cell_type

print (atom_types) 

atm2GPa = 0.000101325
A2Bohr = 1.889725989
kcalmol_2_Ha = 0.00159362
kcalmol_A_2_Ha_Bohr = kcalmol_2_Ha / A2Bohr

f  = open(file_dump ,"rt")
nconf = 0

fname = "lammps.xyz"
f3 = open( fname, "w")

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
        #print ("NONORTHO")

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

    tmp = f.readline().split()
    tmp = np.array(tmp)
    try:
        idx  = np.where(tmp=='x')[0][0]
        idy  = np.where(tmp=='y')[0][0]
        idz  = np.where(tmp=='z')[0][0]
    except:
        idx  = np.where(tmp=='xs')[0][0]
        idy  = np.where(tmp=='ys')[0][0]
        idz  = np.where(tmp=='zs')[0][0]
    idi  = np.where(tmp=='type')[0][0]

    #print (idq, idfx, idfy, idfz, idx, idy, idz, idi)
    x  = np.zeros(shape=(natom))
    y  = np.zeros(shape=(natom))
    z  = np.zeros(shape=(natom))
    strs = ["aaaa" for x in range(natom)]
    atype = np.array(strs)
    #atype = np.empty(natom,dtype=str)
    for i in range(natom):
        tmp = f.readline().split()
        idt = int(tmp[0])
        x[idt-1]  = float(tmp[idx-2])
        y[idt-1]  = float(tmp[idy-2])
        z[idt-1]  = float(tmp[idz-2])
        iid = int(tmp[idi-2])
        atype[idt-1] = atom_types[iid-1]
        #print (iid, atom_types[iid-1], atype[idt-1])

    f3.write("%d\n" %(natom))
    if (cell_type=="NON_ORTHO"):
        f3.write("%s" %("NON_ORTHO" ))
        f3.write("%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f" %( lx,0,0, xy,ly,0, xz,yz,lz ))
        f3.write("\n")
    elif (cell_type=="cell_9"):
        f3.write("%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f" %( lx,0,0, xy,ly,0, xz,yz,lz ))
        f3.write("\n")
    elif (cell_type=="cell_3"):
        f3.write("%15.6f %15.6f %15.6f" %( lx,ly,lz ))
        f3.write("\n")
    else:
        print ("unknown option")
        exit()

    for i in range(natom):
        f3.write("%s" %(atype[i]))
        f3.write("%15.9f" %(x[i]))
        f3.write("%15.9f" %(y[i]))
        f3.write("%15.9f" %(z[i]))
        f3.write("\n")

    nconf += 1

f.close()
f3.close()

print ("number of frame ", nconf)
