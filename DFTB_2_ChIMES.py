import numpy as np

import argparse
parser = argparse.ArgumentParser(description='DFTB_2_xyzf')
# Arguments supported by the code.
parser.add_argument("--file_dftb_gen",                         default='geo_end.gen',  help='file DFTB gen')
parser.add_argument("--file_dftb_detail",                      default='detailed.out', help='file DFTB detailed')
parser.add_argument("--cell_type",                             default='cell_3',       help='cell_3/cell_9/NON_ORTHO')
parser.add_argument("--export_stress",    action='store_true',                         help='export stress')


args = parser.parse_args()
file_dftb_gen    = args.file_dftb_gen
file_dftb_detail = args.file_dftb_detail
cell_type        = args.cell_type
export_stress    = args.export_stress


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
#print (AtomList)
#print (xyz)

def read_dftb_detail(file_dftb_detail, natom):

    f = open(file_dftb_detail, 'rt')

    while True:

        line = f.readline()

        if line == '': break

        keywords = "Total Forces"
        if keywords in line:
            print ("*****")
            print ("read atomic forces")
            print (line)
            fxyz = np.array([])
            for i in range(natom):
                line = f.readline().split()
                tmp = [float(line[i]) for i in range(1,4)]
                fxyz = np.append(fxyz, np.array(tmp))
                # Unit of forces Ha/Bohr, ChIMES uses this unit

        keywords = "Force related energy:"
        if keywords in line:
            print ("*****")
            print ("read energy")
            print (line)
            eV2kcalmol =  23.0609
            energy = float( line.split()[-2] ) * eV2kcalmol

        keywords = "Total stress tensor"
        if keywords in line:
            print ("*****")
            print ("read stresses")
            print (line)
            stress = np.array([])
            for i in range(3):
                line = f.readline().split()
                tmp = [float(line[i]) for i in range(3)]
                stress = np.append(stress, np.array(tmp))
            au2GPa = 1.0/(0.339893208050290 * 0.0001)
            stress = stress * au2GPa

    fxyz = fxyz.reshape((natom,3))
    stress = stress.reshape((3,3))

    f.close()
    return fxyz, energy, stress

fxyz, energy, stress = read_dftb_detail(file_dftb_detail, natom)

def write_ChIMES_format( natom, AtomList, xyz, cell_xyz, cell_type, fxyz, energy, stress, export_stress):
    f2 = open("DFTB.xyzf", "w")
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
    # stress
    if export_stress:
        for i in range(3):
            f2.write("%15.9f" %( stress[i,i]))
        #xy , xy, yz
        f2.write("%15.9f" %( stress[0,1]))
        f2.write("%15.9f" %( stress[0,2]))
        f2.write("%15.9f" %( stress[1,2]))
    # energy
    f2.write("%20.9f" %( energy ))
    f2.write("\n")

    for i in range(natom):
        f2.write("%s" %( AtomList[i] ))
        for j in range(3):
            f2.write("%15.9f" %( xyz[i,j]))
        for j in range(3):
            f2.write("%15.9f" %( fxyz[i,j]))
        f2.write("\n")
    
    f2.close()

write_ChIMES_format( natom, AtomList, xyz, cell_xyz, cell_type, fxyz, energy, stress, export_stress)
