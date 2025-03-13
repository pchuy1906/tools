import numpy as np
from ase.io import read, write
from ase.geometry import  get_distances
#
#
import argparse
parser = argparse.ArgumentParser(description='Identify the molecule')
# Arguments supported by the code.
parser.add_argument("--file_POSCAR",                  default='POSCAR',   help='file_input vasp (POSCAR)')
parser.add_argument("--rcut",                         default=1.6,        help='rcut(1.5)')

args        = parser.parse_args()
file_POSCAR        = args.file_POSCAR
rcut               = args.rcut

atoms = read(file_POSCAR)
cell = atoms.get_cell()

natom = len(atoms)

tmp_arr = []
for i in range(natom):
    TMP = get_distances(atoms.positions[i], atoms.positions, pbc=True, cell=cell)
    tmp_arr.append(TMP[1][0])
    #print (TMP[1][0])

dist_2d_arr = np.stack(tmp_arr)
dist_2d_arr = (dist_2d_arr+dist_2d_arr.T)/2

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
condensed_distance_matrix = squareform(dist_2d_arr)
linkage_matrix = linkage(condensed_distance_matrix, method='single')
clusters = fcluster(linkage_matrix, rcut, criterion='distance')
print(clusters)
nclusters = np.max(clusters)
print (np.sort(clusters))

print ("number of clusters is ", nclusters)

element_symbols = atoms.get_chemical_symbols()
element_symbols = np.array(element_symbols)
print(element_symbols)
print(type(element_symbols))


def write_xyz(fname, natomi, tmp_cell, atomList, xyz):
    f2 = open( fname, "w")
    f2.write("%1d\n" %( natomi ))
    for i in range(3):
        for j in range(3):
            f2.write("%15.9f " %( tmp_cell[i][j] ))
    f2.write("\n")
    for i in range(natomi):
        f2.write("%s " %( atomList[i] ))
        for j in range(3):
            f2.write("%15.9f " %( xyz[i][j] ))
        f2.write("\n")

tmp_cell = cell[:]
for i in range(nclusters):
    iloc = np.where(clusters == i+1)[0]
    fname = "molecule-"+str(i+1)+".xyz"
    natomi = len(iloc)
    atomList = element_symbols[iloc]
    xyz = atoms.positions[iloc]
    write_xyz(fname, natomi, tmp_cell, atomList, xyz)





