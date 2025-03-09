import numpy as np

import argparse
parser = argparse.ArgumentParser(description='VASP_2_ChIMES')
# Arguments supported by the code.
parser.add_argument("--file_VASP_OUTCAR", default='OUTCAR', help='file OUTCAR')
parser.add_argument("--file_VASP_POSCAR", default='POSCAR', help='file POSCAR')
parser.add_argument("--cell_type",        default='cell_3', help='cell_3/NON_ORTHO')
parser.add_argument("--export_stress",    default=False, action="store_true", help="export stress")

args        = parser.parse_args()
file_VASP_OUTCAR       = args.file_VASP_OUTCAR
file_VASP_POSCAR       = args.file_VASP_POSCAR
cell_type              = args.cell_type
export_stress          = args.export_stress

print ()
print ("Program running")
print ()

# parameters
eV2kcalmol=23.061
eV2Ha=0.0367493
Ang2Bohr=1.889725989
kB2GPa=0.1

print ("***************************************")
print ("reading file POSCAR:", file_VASP_POSCAR)
f = open(file_VASP_POSCAR, 'rt')
readPOSCAR = True
while readPOSCAR:
    for i in range(5):
        line = f.readline()
    line = f.readline().split()
    asym = line
    ntype= len(asym)
    print (line)
    line = f.readline().split()
    anum = [int(line[i]) for i in range(ntype)]
    print (line)
    readPOSCAR = False
f.close()
#print (anum)

AtomList = []
for i in range(ntype):
    ntypei = anum[i]
    for j in range(ntypei):
        AtomList.append(asym[i])

natom = len(AtomList)
xyz = np.zeros(shape=(natom,3))
fxyz = np.zeros(shape=(natom,3))

print ("***************************************")
print ("reading file OUTCAR:", file_VASP_OUTCAR)
f = open(file_VASP_OUTCAR, 'rt')
while True:
    line = f.readline()
    if line == '': break

    keywords = "in kB"
    if keywords in line:
        print ("*****")
        print ("read stresses")
        line = line.split()
        #print (line)
        # VASP order is: XX YY ZZ XY YZ ZX
        stresses = [float(line[i])*kB2GPa for i in range(2,8)]
        # ChIMES input wants: XX YY ZZ XY ZX YZ
        tmp = stresses[4]
        stresses[4] = stresses[5]
        stresses[5] = tmp
        str_stresses = ' '.join(map(str, stresses))
        print (str_stresses)
        done_read_stress = True

    keywords = "direct lattice vectors"
    if keywords in line:
        cell_xyz = np.array([])
        print ("*****")
        print ("read cell parameter")
        for i in range(3):
            line = f.readline().split()
            #print (line)
            tmp = [float(line[i]) for i in range(0,3)]
            cell_xyz = np.append(cell_xyz, np.array(tmp))
        str_cell_xyz = ' '.join(map(str, cell_xyz))
        print (str_cell_xyz)

    keywords = "POSITION"
    if keywords in line:
        print ("*****")
        print ("read atomic coordinates and forces")
        line = f.readline().split()
        for i in range(natom):
            line = f.readline().split()
            #print (line)
            xyz[i,:] = [float(line[i]) for i in range(0,3)]
            fxyz[i,:] = [float(line[i]) for i in range(3,6)]
            ''' VASP has force (eV/A), 
              while ChIMES used (Ha/Bohr) for input
              and (kcal/mol/A) for output
            '''
        fxyz = fxyz * eV2Ha/(Ang2Bohr)

    keywords = "free  energy"
    if keywords in line:
        print ("*****")
        print ("read energy")
        #print (line)
        energy = float(line.split()[4])
        print (energy)
        energy *= eV2kcalmol
f.close()

f2 = open("VASP_2_ChIMES.xyzf", "w")

f2.write("%1d\n" %( natom ))

# write the cell parameter
if (cell_type=="cell_3"):
    for i in range(3):
        j = 4*i
        f2.write("%15.9f" %( cell_xyz[j]))
elif (cell_type=="NON_ORTHO"):
    f2.write("%s" %( "NON_ORTHO " ))
    for i in range(9):
        f2.write("%15.9f" %( cell_xyz[i]))
# write stresses
if export_stress:
    for i in range(6):
        f2.write("%15.9f" %( stresses[i]))
# write energy
f2.write("%20.9f" %( energy ))

f2.write("\n")
for i in range(natom):
    f2.write("%s" %( AtomList[i] ))
    for j in range(3):
        f2.write("%15.9f" %( xyz[i,j]))
    for j in range(3):
        f2.write("%15.9f" %( fxyz[i,j]))
    f2.write("\n")

