import sys
module_path = "/g/g92/pham20/tools/others"
sys.path.append(module_path)
import utils

import numpy as np
 
import argparse
parser = argparse.ArgumentParser(description='convert cell parameter to lengths and angles')
# Arguments supported by the code.
parser.add_argument("--file_cell_xyz",      default='cell.dat',      help='cellXYZ/data.lammps/cell.dat')
parser.add_argument("--cell_type",          default='cell_3',        help='cellXYZ/lammps_data/cell_9')
args        = parser.parse_args()
file_cell_xyz = args.file_cell_xyz
cell_type     = args.cell_type


print ("")
print ("read cell XYZ")
if (cell_type=="cell_9"):
    cellxyz = np.loadtxt(file_cell_xyz)
    cellxyz = cellxyz.reshape(3, 3)
    print (cellxyz)
elif (cell_type=="lammps_data"):
    f = open(file_cell_xyz, 'rt')
    while True:
        line = f.readline()
        if line == '': break
        keywords = "xlo xhi"
        if keywords in line:
            lo_hi = line.split()[:2]
            x = float(lo_hi[1]) - float(lo_hi[0])
        keywords = "ylo yhi"
        if keywords in line:
            lo_hi = line.split()[:2]
            y = float(lo_hi[1]) - float(lo_hi[0])
        keywords = "zlo zhi"
        if keywords in line:
            lo_hi = line.split()[:2]
            z = float(lo_hi[1]) - float(lo_hi[0])
        keywords = "xy xz yz"
        if keywords in line:
            tmp = line.split()[:3]
            xy,xz,yz = [float(i) for i in tmp]

    cellxyz = np.array([[x, 0, 0],
                        [xy, y, 0],
                        [xz, yz, z]])
            


Volume = np.linalg.det(cellxyz)
print ("unit-cell volume:",Volume)

print ("")
print ("convert cellxyz to lengths and angles")
lengths, angles = utils.Cell_XYZ_ABC(cellxyz)
#print (lengths,angles)
formatted = [f"{x:.2f}" for x in lengths] + [f"{x:.2f}" for x in angles]
print(' '.join(formatted))
