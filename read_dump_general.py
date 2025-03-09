import numpy as np
 
# read input 
import argparse
parser = argparse.ArgumentParser(description='split file xyz')
# Arguments supported by the code.
parser.add_argument("--file_dump",   default='file.dump', help='file with format dump lammps')
parser.add_argument("--file_energy", default='file.ener', help='file with format log lammps')
parser.add_argument("--export_quan", default='xyzfes',    help='export xyzfe/xyzfes')
parser.add_argument("--cell_type",   default='NON_ORTHO', help='cell_type (cell_3/cell_9/NON_ORTHO)')
parser.add_argument('--atom_types',  nargs='+')
parser.add_argument("--lmp_unit",    default='real',    help='real/metal')
parser.add_argument("--read_charge", default=False, action="store_true", help="read point charge")

args    = parser.parse_args()
file_dump     = args.file_dump
file_energy   = args.file_energy
atom_types    = args.atom_types
export_quan   = args.export_quan
cell_type     = args.cell_type
lmp_unit      = args.lmp_unit
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

if (export_quan=="xyzfes"):
    fname = "LAMMPS.xyzfes"
elif (export_quan=="xyzfe"):
    fname = "LAMMPS.xyzfe"
else:
    print ("unknown option")
    exit()
f3 = open( fname, "w")

f0  = open(file_energy ,"rt")
tmp = f0.readline().split()
tmp = np.array(tmp)
ipe  = np.where(tmp=='PotEng')[0][0]
ixx  = np.where(tmp=='Pxx')[0][0]
iyy  = np.where(tmp=='Pyy')[0][0]
izz  = np.where(tmp=='Pzz')[0][0]
ixy  = np.where(tmp=='Pxy')[0][0]
ixz  = np.where(tmp=='Pxz')[0][0]
iyz  = np.where(tmp=='Pyz')[0][0]

while True:

    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break
    tmp  = f.readline()
    tmp  = f.readline()
    tmp  = f.readline()
    natom = int(tmp)

    tmp = f.readline()
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

    tmp = f.readline().split()
    tmp = np.array(tmp)
    if (read_charge):
        idq  = np.where(tmp=='q')[0][0]
    idfx = np.where(tmp=='fx')[0][0]
    idfy = np.where(tmp=='fy')[0][0]
    idfz = np.where(tmp=='fz')[0][0]
    idx  = np.where(tmp=='x')[0][0]
    idy  = np.where(tmp=='y')[0][0]
    idz  = np.where(tmp=='z')[0][0]
    idi  = np.where(tmp=='type')[0][0]

    #print (idq, idfx, idfy, idfz, idx, idy, idz, idi)
    q  = np.zeros(shape=(natom))
    x  = np.zeros(shape=(natom))
    y  = np.zeros(shape=(natom))
    z  = np.zeros(shape=(natom))
    fx = np.zeros(shape=(natom))
    fy = np.zeros(shape=(natom))
    fz = np.zeros(shape=(natom))
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
        if (lmp_unit == "real"):
            f_convert = kcalmol_A_2_Ha_Bohr
            p_convert = atm_2_GPa
            e_convert = 1.0
        elif (lmp_unit == "metal"):
            f_convert = eV_A_2_Ha_Bohr
            p_convert = bar_2_GPa
            e_convert = eV_2_kcalmol


        fx[idt-1] = float(tmp[idfx-2]) * f_convert
        fy[idt-1] = float(tmp[idfy-2]) * f_convert
        fz[idt-1] = float(tmp[idfz-2]) * f_convert
        iid = int(tmp[idi-2])
        atype[idt-1] = atom_types[iid-1]

    TPM  = f0.readline().split()
    energy = float(TPM[ipe]) * e_convert

    Pxx = float(TPM[ixx]) * p_convert 
    Pyy = float(TPM[iyy]) * p_convert 
    Pzz = float(TPM[izz]) * p_convert 
    Pxy = float(TPM[ixy]) * p_convert 
    Pxz = float(TPM[ixz]) * p_convert 
    Pyz = float(TPM[iyz]) * p_convert 

    f3.write("%d\n" %(natom))
    if (export_quan=="xyzfes"):
        if (cell_type=="NON_ORTHO"):
            f3.write("%s" %("NON_ORTHO" ))
            f3.write("%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f" %( lx,0,0, 0,ly,0, 0,0,lz ))
            f3.write("%20.6f %20.6f %20.6f %20.6f %20.6f %20.6f" %( Pxx, Pyy, Pzz, Pxy, Pxz, Pyz ))
            f3.write("%20.6f" %( energy ))
            f3.write("\n")
        elif (cell_type=="cell_3"):
            f3.write("%15.6f %15.6f %15.6f" %( lx,ly,lz ))
            f3.write("%20.6f %20.6f %20.6f %20.6f %20.6f %20.6f" %( Pxx, Pyy, Pzz, Pxy, Pxz, Pyz ))
            f3.write("%20.6f" %( energy ))
            f3.write("\n")
        else:
            print ("unknown option")
            exit()
    elif (export_quan=="xyzfe"):
        f3.write("%15.6f %15.6f %15.6f" %( lx,ly,lz ))
        f3.write("%20.6f" %( energy ))
        f3.write("\n")
    else:
        print ("unknown option")
        exit()


    for i in range(natom):
        if (read_charge):
            f2.write("%15.9f\n" %(q[i]))
        f3.write("%s" %(atype[i]))
        f3.write("%15.9f" %(x[i]))
        f3.write("%15.9f" %(y[i]))
        f3.write("%15.9f" %(z[i]))
        f3.write("%15.9f" %(fx[i]))
        f3.write("%15.9f" %(fy[i]))
        f3.write("%15.9f" %(fz[i]))
        f3.write("\n")


    nconf += 1


f2.close()


