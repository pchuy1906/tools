import numpy as np

import argparse
parser = argparse.ArgumentParser(description='gen_2_xyz')
# Arguments supported by the code.
parser.add_argument("--file_dftb_gen",                         default='geo_end.gen',  help='file DFTB gen')
parser.add_argument("--cell_type",                             default='cell_3',       help='cell_3/cell_9/NON_ORTHO')

args = parser.parse_args()
file_dftb_gen    = args.file_dftb_gen
cell_type        = args.cell_type

def read_dftb_gen(file_dftb_gen):

    f = open(file_dftb_gen, 'rt')

    line = f.readline().split()
    natom = int(line[0])
    type_xyz = line[1]
    print ("number of atom = %d" %natom)

    atypes = f.readline().split()
    xyz = np.array([])
    cell_xyz = np.array([])
    AtomList = []

    for i in range(natom):

        line = f.readline().split()

        AtomList.append( atypes[ int(line[1])-1 ] )

        tmp = [float(line[j]) for j in range(2,5)]
        xyz = np.append(xyz, np.array(tmp))

    line = f.readline()

    line = f.readline().split()
    tmp = [float(line[j]) for j in range(3)]
    cell_xyz = np.append(cell_xyz, np.array(tmp))

    line = f.readline().split()
    tmp = [float(line[j]) for j in range(3)]
    cell_xyz = np.append(cell_xyz, np.array(tmp))

    line = f.readline().split()
    tmp = [float(line[j]) for j in range(3)]
    cell_xyz = np.append(cell_xyz, np.array(tmp))

    f.close()

    cell_xyz = cell_xyz.reshape((3,3))
    xyz = xyz.reshape((natom,3))
    if (type_xyz=="F"):
        xyz = np.dot(xyz,cell_xyz)

    return natom, AtomList, xyz, cell_xyz

natom, AtomList, xyz, cell_xyz = read_dftb_gen(file_dftb_gen)

def write_xyz_format( natom, AtomList, xyz, cell_xyz, cell_type):
    f2 = open("gen_2_xyz.xyz", "w")
    f2.write("%1d\n" %( natom ))
    # cell parameters
    if (cell_type == "cell_3"):
        for i in range(3):
            f2.write("%15.9f" %( cell_xyz[i,i]))
    elif (cell_type == "cell_9"):
        for i in range(3):
            for j in range(3):
                f2.write("%15.9f" %( cell_xyz[i,j]))
    elif (cell_type == "NON_ORTHO"):
        f2.write("NON_ORTHO ")
        for i in range(3):
            for j in range(3):
                f2.write("%15.9f" %( cell_xyz[i,j]))
    else:
        print ("unknown cell_type")
    f2.write("\n")

    for i in range(natom):
        f2.write("%s" %( AtomList[i] ))
        for j in range(3):
            f2.write("%15.9f" %( xyz[i,j]))
        f2.write("\n")
    
    f2.close()

write_xyz_format( natom, AtomList, xyz, cell_xyz, cell_type)
