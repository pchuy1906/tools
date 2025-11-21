import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np

import argparse
parser = argparse.ArgumentParser(description='QE_2_xyzf')
# Arguments supported by the code.
parser.add_argument("--file_QE_output",         default='output', help='file QE output')
parser.add_argument("--cell_type",              default='cell_3', help='cell_3/NON_ORTHO')
parser.add_argument("--export_stress",          default=False, action="store_true", help="export stresses")


args        = parser.parse_args()
file_QE_output        = args.file_QE_output
cell_type             = args.cell_type
export_stress         = args.export_stress

print ("reading file_QE_output:", file_QE_output)

xyz = np.array([])
cell_xyz = np.array([])
fxyz = np.array([])
stresses = np.array([])
AtomList = []

f = open(file_QE_output, 'rt')
while True:
    line = f.readline()
    if line == '': break

    keywords = "lattice parameter (alat)"
    if keywords in line:
        print (line)
        alat = float(line.split()[4])

    keywords = "number of atoms/cell"
    if keywords in line:
        print ("*****")
        print ("read the number of atom")
        print (line)
        natom = line.split()[-1]
        natom = int(natom)
        print (natom)

    keywords = "crystal axes:"
    if keywords in line:
        print ("*****")
        print ("read cell parameter")
        for i in range(3):
            line = f.readline().split()
            #print (line)
            tmp = [float(line[i]) for i in range(3,6)]
            cell_xyz = np.append(cell_xyz, np.array(tmp))

    keywords = "positions (alat units)"
    if keywords in line:
        print ("*****")
        print ("read atomic coordinates")
        print (line)
        for i in range(natom):
            line = f.readline().split()
            #print (line)
            tmp = [float(line[i]) for i in range(6,9)]
            #print (tmp)
            xyz = np.append(xyz, np.array(tmp))
            AtomList.append(line[1])

    keywords = "!    total energy"
    if keywords in line:
        print ("*****")
        print ("read energy")
        print (line)
        energy = float(line.split()[4])

    keywords = "Forces acting on atoms"
    if keywords in line:
        print ("*****")
        print ("read atomic forces")
        print (line)
        line = f.readline()
        for i in range(natom):
            line = f.readline().split()
            #print (line)
            tmp = [float(line[i]) for i in range(6,9)]
            fxyz = np.append(fxyz, np.array(tmp))

    keywords = "kbar"
    if keywords in line:
        print ("*****")
        print ("read stresses")
        print (line)
        for i in range(3):
            line = f.readline().split()
            #print (line)
            tmp = [float(line[i]) for i in range(3,6)]
            stresses = np.append(stresses, np.array(tmp))

f.close()

#print (cell_xyz)
cell_xyz = cell_xyz*alat/1.889726
cell_9 = cell_xyz
cell_xyz = cell_xyz.reshape((3,3))
#print (cell_xyz)

xyz = xyz.reshape((natom,3))
xyz = xyz*alat/1.889726

fxyz = fxyz.reshape((natom,3))
#QE output (Ry/Bohr)
#ChIMES input (Ha/Bohr)
Ry2Ha = 0.5
fxyz = fxyz * Ry2Ha

Ry2kcalmol = 313.75470835207074 
energy = energy * Ry2kcalmol

stresses = stresses * 0.1
stresses = stresses.reshape((3,3))

fname = "QE.xyzf"
utils.write_xyzf(fname, natom, AtomList, xyz, cell_xyz, cell_type, fxyz, energy, stresses, export_stress)

