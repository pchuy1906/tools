import numpy as np

import argparse
parser = argparse.ArgumentParser(description='export FORCES, ENERGIES, and STRESS of file xyzfes')
# Arguments supported by the code.
parser.add_argument("--file_xyz",           default='file.xyz',                help='file_XYZ format xyz')
parser.add_argument("--file_atomic_energy", default='atomic_energies.dat',     help='input/output the atomic energies')
parser.add_argument("--fit_atomic_energy", action='store_true',                help='fit or use atomic energies')
parser.add_argument("--cell_type",          default='',                        help='cell_3/cell_9/NON_ORTHO')
parser.add_argument('--atom_type', nargs='+')

args        = parser.parse_args()
file_xyz           = args.file_xyz
file_atomic_energy = args.file_atomic_energy
fit_atomic_energy  = args.fit_atomic_energy
cell_type          = args.cell_type
atom_type          = args.atom_type

print ("-----------------------------------------------------------")
print ("")
print ("")
print ("read fileXYZ:", file_xyz)
f  = open(file_xyz ,"rt")

def number_CHNO(atom_type, molecule_order, molecule_natom):
    atom_type = np.array(atom_type)
    res = np.zeros(len(atom_type))
    for id_tmp in range(len(molecule_order)):
        atom_loc = np.where(atom_type == molecule_order[id_tmp])
        res[atom_loc] = molecule_natom[id_tmp]
    return res

istruc = 0

Amatrix = np.array([])
bmatrix = []

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
    bmatrix.append(float(tmp[ncomment-1]))

    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(0,natom):
        tmp = f.readline()
        tmp = tmp.split()
        atomList.append(tmp[0])
        xyz[k,0],xyz[k,1],xyz[k,2] =  float(tmp[1]),float(tmp[2]),float(tmp[3])
        force[3*k+0],force[3*k+1],force[3*k+2] = float(tmp[4]),float(tmp[5]),float(tmp[6])
    molecule_order, molecule_natom = np.unique(atomList, return_counts=True)
    Amatrix = np.append(Amatrix, number_CHNO(atom_type, molecule_order, molecule_natom) )
    istruc += 1
f.close

Amatrix = Amatrix.reshape((istruc, len(atom_type)))
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


f2 = open("subtract_atomic_energy_"+file_xyz, "w")

f  = open(file_xyz ,"rt")
istruc=0

def write_comment_line(f2, cell_type, cell_9, stress, energy):
    if cell_type=="NON_ORTHO":
        f2.write("%s " %("NON_ORTHO"))
        for i in range(9):
            f2.write("%12.4f " %(cell_9[i]))
        if (len(stress)==6):
            for i in range(6):
                f2.write("%12.4f " %(stress[i]))
        f2.write("%20.6f " %(energy))
    elif cell_type=="cell_3":
        for i in [0, 4, 8]:
            f2.write("%12.4f " %(cell_9[i]))
        if (len(stress)==6):
            for i in range(6):
                f2.write("%12.4f " %(stress[i]))
        f2.write("%20.6f " %(energy))
    else:
        print ("unknown cell_type", cell_type)
    f2.write("\n" )

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp)
    f2.write(tmp)

    tmp = f.readline().split()

    ncomment = len(tmp)
    energy = (float(tmp[ncomment-1]) - np.dot(Amatrix[istruc,:],Atomic_energies[:]))

    if (ncomment==17):
        # format: NON_ORTHO cell(1,1:3), cell(2,1:3), cell(3,1:3), stress(1:6), energy
        cell_9 = [float(tmp[i]) for i in range(1,10)]
        stress = [float(tmp[i]) for i in range(10,16)]
    elif (ncomment==11):
        # format: NON_ORTHO cell(1,1:3), cell(2,1:3), cell(3,1:3), energy
        cell_9 = [float(tmp[i]) for i in range(1,10)]
        stress = []
    elif (ncomment==10):
        # format: cell(1), cell(2), cell(3), stress(1:6), energy
        cell_9 = [0.0 for i in range(1,10)]
        cell_9[0] = float(tmp[0])
        cell_9[4] = float(tmp[1])
        cell_9[8] = float(tmp[2])
        stress = [float(tmp[i]) for i in range(3,9)]
    elif (ncomment==4):
        # format: cell(1), cell(2), cell(3), energy
        cell_9 = [0.0 for i in range(1,10)]
        cell_9[0] = float(tmp[0])
        cell_9[4] = float(tmp[1])
        cell_9[8] = float(tmp[2])
        stress = []
    else:
        print ("unknown option")
        exit()

    write_comment_line(f2, cell_type, cell_9, stress, energy)

    for k in range(natom):
        tmp = f.readline()
        f2.write(tmp)
    istruc+=1
f.close

