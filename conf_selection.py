import numpy as np

Ha2kcalmol = 627.509
Bohr2Angstrom = 0.529177

import argparse
parser = argparse.ArgumentParser(description='group same molecules to one xyz file')
# Arguments supported by the code.
parser.add_argument("--file_xyz",                 default='file.xyz',     help='file_XYZ format xyz')
parser.add_argument("--file_out_xyz",             default='file_out.xyz', help='file_out_XYZ format xyz')
parser.add_argument("--icell",        type=int,   default=1,              help='1/cell_3 2/cell_9 3/NON_ORTHO')
parser.add_argument("--fcut",         type=float, default=2000.0,         help='maximum value of force to be considered')


args        = parser.parse_args()
file_xyz           = args.file_xyz
file_out_xyz       = args.file_out_xyz
icell              = args.icell
fcut               = args.fcut

print ("-----------------------------------------------------------")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"rt")

def number_CHNO(atom_order, molecule_order, molecule_natom):
    res = np.array([0,0,0,0])
    for id_tmp in range(len(molecule_order)):
        atom_loc = np.where(atom_order == molecule_order[id_tmp])
        res[atom_loc] = molecule_natom[id_tmp]
    return res

def iname(molecule_order, molecule_natom):
    nlen = len(molecule_order)
    molecule_name = ""
    for i in range(nlen):
        aaa = molecule_order[i] + str(molecule_natom[i]) + "_"
        molecule_name = molecule_name + aaa
    return molecule_name


istruc = 0
atom_order = np.array(['C','H','N','O'])

f2  = open(file_out_xyz ,"w")

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break
    #v, e = [float(x) for x in line.split()[:2]]

    natom = int(tmp)
    
    tmp = f.readline().split()
    if (icell==3):
        cell_9 = [float(x) for x in tmp[1:10]]
        stress_6 = [float(x) for x in tmp[10:16]]
        energy = tmp[16]
    elif (icell==1):
        cell_3 = [float(x) for x in tmp[0:3]]
        energy = tmp[3]
    else:
        print ("unknown option !! error !!")


    atomList = []
    xyz = np.zeros(shape=(natom,3))
    fxyz = np.zeros(shape=(natom,3))

    iwrite = True
    for k in range(natom):
        tmp = f.readline().split()
        atomList.append(tmp[0])
        xyz[k,0], xyz[k,1], xyz[k,2] = float(tmp[1]),float(tmp[2]),float(tmp[3])
        fxyz[k,0],fxyz[k,1],fxyz[k,2] = float(tmp[4]),float(tmp[5]),float(tmp[6])

    #print (np.amax(np.abs(fxyz)))
    f_max = np.amax(np.abs(fxyz)) * Ha2kcalmol/Bohr2Angstrom
    #print (f_max)
    molecule_order, molecule_natom = np.unique(atomList, return_counts=True)
    molecule_name = iname(molecule_order, molecule_natom)
    #print (molecule_name)

    istruc += 1
    if (f_max < fcut):
        f2.write("%1d\n" %( natom ))
        if (icell==3):
            #f2.write("%s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %s\n" %( "NON_ORTHO", cell_9[0], cell_9[1], cell_9[2], cell_9[3], cell_9[4], cell_9[5], cell_9[6], cell_9[7], cell_9[8], stress_6[0], stress_6[1], stress_6[2], stress_6[3], stress_6[4], stress_6[5], energy))

            f2.write("%s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f \
                     %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %s\n" \
                     %( "NON_ORTHO", cell_9[0], cell_9[1], cell_9[2], cell_9[3], cell_9[4], cell_9[5], cell_9[6], cell_9[7],cell_9[8], \
                     stress_6[0], stress_6[1], stress_6[2], stress_6[3], stress_6[4], stress_6[5], energy))
        elif (icell==1):
            f2.write("%15.9f %15.9f %15.9f %s\n" %( cell_3[0], cell_3[1], cell_3[2], energy) )
        else:
            print ("unknown option !! error !!")


        for k in range(natom):
            f2.write("%s %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n" %(atomList[k], xyz[k,0], xyz[k,1], xyz[k,2], fxyz[k,0],fxyz[k,1],fxyz[k,2]))


f.close

