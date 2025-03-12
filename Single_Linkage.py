import numpy as np
from ase.io import read, write
from ase.geometry import  get_distances


import argparse
parser = argparse.ArgumentParser(description='Identify the molecule')
# Arguments supported by the code.
parser.add_argument("--file_POSCAR",                  default='POSCAR',   help='file_input vasp (POSCAR)')
parser.add_argument("--rcut",                         default=1.6,        help='rcut(1.5)')

args        = parser.parse_args()
file_POSCAR        = args.file_POSCAR
rcut               = args.rcut


atoms = read(file_POSCAR)
natom = len(atoms)

tmp_arr = []
for i in range(natom):
    TMP = get_distances(atoms.positions[i], atoms.positions)
    tmp_arr.append(TMP[1][0])
    #print (TMP[1][0])

dist_2d_arr = np.stack(tmp_arr)
print (dist_2d_arr)


#from sklearn.cluster import AgglomerativeClustering
#clusterer = AgglomerativeClustering(n_clusters=None, metric="precomputed", linkage="average", distance_threshold=1.5)
#clusters = clusterer.fit_predict(dist_2d_arr)
#print(clusters)

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


condensed_distance_matrix = squareform(dist_2d_arr)
print (condensed_distance_matrix)
linkage_matrix = linkage(condensed_distance_matrix, method='single')
clusters = fcluster(linkage_matrix, rcut, criterion='distance')
print(clusters)
print (np.sort(clusters))



#        Z = linkage(dist_4_linkage, 'single')
#        clusters = fcluster(Z, rcut, criterion='distance')

