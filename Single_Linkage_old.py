import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
#from matplotlib import pyplot as plt
from sys import exit

import argparse
parser = argparse.ArgumentParser(description='Identify the molecule')
# Arguments supported by the code.
parser.add_argument("--file_xyz",                     default='file.xyz', help='file_input (file.xyz)')
parser.add_argument("--rcut",                         default=1.5,        help='rcut(1.5)')
parser.add_argument("--natom_per_molecule", type=int, default=0,          help='number of atom per molecule, if specify, skip the cluster')
parser.add_argument("--print_molecules",    type=int, default=0,          help='print all molecules (0/1)')
parser.add_argument("--cell_option",                  default='cell_3',   help='cell_option = "NON_ORTHO", "cell_3", or "cell_9" ')
parser.add_argument("--file_cluster",                 default='abcdef',   help='file cluster')

args        = parser.parse_args()
file_xyz           = args.file_xyz
rcut               = args.rcut
natom_per_molecule = args.natom_per_molecule
print_molecules    = args.print_molecules
cell_option        = args.cell_option
file_cluster       = args.file_cluster

f  = open(file_xyz ,"r")
natom = f.readline()
natom = int(natom)
print ("The number of atoms is: %d" % natom)

tmp = f.readline().split()
cell_xyz = np.zeros(shape=(3,3))
if (cell_option == "cell_3"):
    for i in range(3):
        cell_xyz[i,i] = float(tmp[i])
elif (cell_option == "cell_9"):
    tmp_cell = [float(tmp[i]) for i in range(9)]
    tmp_cell = np.array(tmp_cell)
    tmp_cell = tmp_cell.reshape((3, 3))
    cell_xyz = tmp_cell
elif (cell_option == "NON_ORTHO"):
    tmp_cell = [float(tmp[i]) for i in range(1,10)]
    tmp_cell = np.array(tmp_cell)
    tmp_cell = tmp_cell.reshape((3, 3))
    cell_xyz = tmp_cell
else:
    print ("cell paramters can not read!!")
    exit(0)

xyz = np.zeros(shape=(natom,3))
AtomList = []

for k in range(natom):
    tmp = f.readline().split()
    AtomList.append(tmp[0])
    xyz[k,0] =  float(tmp[1])
    xyz[k,1] =  float(tmp[2])
    xyz[k,2] =  float(tmp[3])
f.close


# ## 2. Calculate distance matrix between atoms
print(cell_xyz)
print()

# compute distance using PBC, working only for orthorhombic cell?
def _distance(x0, x1, cell_xyz):
    cell = [cell_xyz[0,0], cell_xyz[1,1], cell_xyz[2,2]]
    cell = np.array(cell)
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * cell, delta - cell, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

def distance(xyz1, xyz2, cell_xyz):
    xyz = xyz1 - xyz2
    xyz = np.dot(xyz, np.linalg.inv(cell_xyz))
    #print (xyz.shape)
    nx = xyz.shape[0]
    ny = xyz.shape[1]
    nxyz = np.zeros(shape=(nx,ny))
    #print (nx, ny, nxyz)
    for i in range(nx):
        for j in range(ny):
            nxyz[i,j] = xyz[i,j]-int(round(xyz[i,j]))
    xyz = np.dot(nxyz, cell_xyz)
    #print (nx, ny, xyz)
    return np.sqrt((xyz ** 2).sum(axis=-1)), xyz1-xyz 


if (natom_per_molecule==0):
    dist = np.array([])
    for i in range(natom):
        tmp_dist, tmp_xyz = distance(xyz[i+1:,:], xyz[i,:], cell_xyz)
        dist = np.append(dist, tmp_dist)
    dist_4_linkage = dist
    print (dist_4_linkage)
    print (len(dist_4_linkage))

# ## 3. Perform single linkage cluster to identify the molecule
if file_cluster == "abcdef":

    if (natom_per_molecule==0):
        Z = linkage(dist_4_linkage, 'single')
        clusters = fcluster(Z, rcut, criterion='distance')
    else:
        clusters = []
        for i in range(natom/natom_per_molecule):
            for j in range(natom_per_molecule):
                clusters.append(i+1)
    clusters = np.array(clusters)
    np.savetxt('cluster.txt', clusters, fmt='%d')

else:

    clusters = np.loadtxt(file_cluster, dtype=int)


nmolecule = clusters.max()
print ('the number of molecules is:', nmolecule)
print (clusters)

def find_center(xyz, index, cell_xyz):
    # find center
    Center = np.zeros(shape=(3))
    for iatom in index:
        for ixyz in range(3):
            Center[ixyz] += xyz[iatom,ixyz]
    Center = Center/len(index)
    # move the center to the cell
    NewCenter = np.dot(Center, np.linalg.inv(cell_xyz))
    for ixyz in range(3):
        NewCenter[ixyz] = NewCenter[ixyz]-np.floor(NewCenter[ixyz])
    NewCenter = np.dot(NewCenter, cell_xyz)
    return Center, NewCenter


# In[7]:


f3 = open('molecule-all.xyz', "w")
f3.write("%10d\n" %( natom ))

for i in range(3):
    for j in range(3):
        f3.write("%15.9f" %( cell_xyz[i,j]))
f3.write("\n")

for imolecule in range(1,nmolecule+1):
    index = [i for i in range(len(clusters)) if clusters[i]==imolecule]
    print (index)
    print ('working with molecule-', imolecule, len(index))
    
    print ('moving all atoms to form a molecule...')
    collect_index = []
    for iatom in range(1,len(index)):
        print ('---work with atom index', index[iatom])
        collect_index.append(index[iatom-1])
        xyz_center,_ = find_center(xyz, collect_index, cell_xyz)
        xyz_center = xyz_center.reshape(1,3)
        print ('mol_center is:', xyz_center)
        print ('old coordinates:', xyz[index[iatom],:])
        _distance, _xyz2 = distance(xyz_center, xyz[index[iatom],:], cell_xyz)
        xyz[index[iatom],:] = _xyz2
        print ('new coordinates:', _xyz2)
    
    print ('calculate molecule center and move the molecule to the cell_xyz')
    Center, NewCenter = find_center(xyz, index, cell_xyz)
    for iatom in index:
        xyz[iatom,:] += NewCenter - Center
    #
    if (print_molecules==1):
        f2 = open('molecule-'+str(imolecule)+'.xyz', "w")
        f2.write("%10d\n" %( len(index) ))
        if cell_option=="cell_9":
            for i in range(3):
                for j in range(3):
                    f2.write("%15.9f" %( cell_xyz[i,j]))
        elif cell_option=="cell_3":
            for i in range(3):
                f2.write("%15.9f" %( cell_xyz[i,i]))
        else:
            print ("wrong cell_option")

        f2.write("\n")
    for k in index:
        if (print_molecules==1):
            f2.write("%4s" %( AtomList[k] ))
        f3.write("%4s" %( AtomList[k] ))
        for ixyz in range(3):
            if (print_molecules==1):
                f2.write("%15.9f" %( xyz[k,ixyz] ))
            f3.write("%15.9f" %( xyz[k,ixyz] ))
        if (print_molecules==1):
            f2.write("\n")
        f3.write("\n")
    if (print_molecules==1):
        f2.close
f3.close

