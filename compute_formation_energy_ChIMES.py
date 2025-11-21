import numpy as np

import argparse
parser = argparse.ArgumentParser(description='compute formation energy')
# Arguments supported by the code.
parser.add_argument("--file_xyz",                         default='file.xyz',   help='file xyz')
parser.add_argument("--file_lammps_log",                  default='log.lammps', help='file log.lammps')
parser.add_argument("--file_AE",                          default='AE.dat',     help='file_AE.dat')

args = parser.parse_args()
file_xyz        = args.file_xyz
file_lammps_log = args.file_lammps_log
file_AE         = args.file_AE


def read_xyz(file_xyz, cell_type):
    f  = open(file_xyz ,"r")
    
    natom = f.readline()
    natom = int(natom)
    print ("the number of atom:", natom)
    
    cell_3_3 = np.zeros(shape=(3,3))
    #print (cell_3_3)
    tmp = f.readline().split()
    
    if cell_type=="cell_3":
        for k in range(3):
            cell_3_3[k,k] = float(tmp[k])
            #print (tmp[k])
            #print (cell_3_3)

    elif cell_type=="cell_9":
        cell_3_3[0,0] = float(tmp[0])
        cell_3_3[0,1] = float(tmp[1])
        cell_3_3[0,2] = float(tmp[2])
    
        cell_3_3[1,0] = float(tmp[3])
        cell_3_3[1,1] = float(tmp[4])
        cell_3_3[1,2] = float(tmp[5])
    
        cell_3_3[2,0] = float(tmp[6])
        cell_3_3[2,1] = float(tmp[7])
        cell_3_3[2,2] = float(tmp[8])
    elif cell_type=="NON_ORTHO":
        cell_3_3[0,0] = float(tmp[1])
        cell_3_3[0,1] = float(tmp[2])
        cell_3_3[0,2] = float(tmp[3])
    
        cell_3_3[1,0] = float(tmp[4])
        cell_3_3[1,1] = float(tmp[5])
        cell_3_3[1,2] = float(tmp[6])
    
        cell_3_3[2,0] = float(tmp[7])
        cell_3_3[2,1] = float(tmp[8])
        cell_3_3[2,2] = float(tmp[9])
    else:
        print ("unknown cell_type", cell_type)
        exit()
    
    #print (cell_3_3)

    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(natom):
        tmp = f.readline()
        tmp = tmp.split()
        atomList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
    f.close
    return natom, cell_3_3, atomList, xyz


def read_lammps_log(file_lammps_log):

    f = open(file_lammps_log, 'rt')

    while True:

        line = f.readline()

        if line == '': break

        keywords = "Energy initial, next-to-last, final ="
        if keywords in line:
            print (line)
            line_next = f.readline().split()
            print (line_next)
            Final_Energy = float(line_next[-1])

    return Final_Energy

print ()
print ()
Etot_LAMMPS = read_lammps_log(file_lammps_log)
print ("Total LAMMPS energy Etot_LAMMPS (kcal/mol) = ", Etot_LAMMPS)
print ()
print ()

print ("read file.xyz")
cell_type = 'cell_3'
natom, cell_3_3, atomList, xyz = read_xyz(file_xyz, cell_type)

def number_CHNO(atom_order, molecule_order, molecule_natom):
    res = np.array([0,0,0,0])
    for id_tmp in range(len(molecule_order)):
        atom_loc = np.where(atom_order == molecule_order[id_tmp])
        res[atom_loc] = molecule_natom[id_tmp]
    return res

atom_order = np.array(['C','H','N','O'])
molecule_order, molecule_natom = np.unique(atomList, return_counts=True)
nCHNO = number_CHNO(atom_order, molecule_order, molecule_natom)
print (nCHNO)

Atomic_energies = np.loadtxt(file_AE)
print ("atomic energies:")
print (Atomic_energies)

print ()
print ()
Etot_VASP = Etot_LAMMPS + np.dot(np.array(nCHNO),np.array(Atomic_energies))
print ("Total VASP energy Etot_VASP (kcal/mol) = ", Etot_VASP)
print ()
print ()

AE_VASP_C = -2422.64418787 / 240
AE_VASP_H = -6.91235731 / 2
AE_VASP_N = -18.57591407 / 2
AE_VASP_O = -10.64265818 / 2

eV2kcalmol = 23.06035
AE_VASP = np.array([AE_VASP_C,AE_VASP_H,AE_VASP_N,AE_VASP_O]) * eV2kcalmol
Eform_VASP = Etot_VASP - np.dot(np.array(nCHNO),np.array(AE_VASP))
print ("formation energy (kcal/mol) = ", Eform_VASP)

