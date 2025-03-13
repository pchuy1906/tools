import numpy as np

# read input 
import argparse
parser = argparse.ArgumentParser(description='compute QEq')
# Arguments supported by the code.
parser.add_argument("--file_Q",             default='charge.dat', help='file charge')
parser.add_argument("--file_xyz",    default='file.xyz', help='file with format xyz/xyzfe/xyzfes')
parser.add_argument("--option_cell", default='cell_3',   help='NON_ORTHO/cell_3/cell_9')

args    = parser.parse_args()
file_Q      = args.file_Q
file_xyz    = args.file_xyz
option_cell = args.option_cell

q = np.loadtxt(file_Q)
#print (q)

def compute_forces_QEq(Q, dist_matrix, xyz, cell_3):
    Coulomb_k = 8.987* 10.0**9.0
    Coulomb_e = 1.60217662* 10.0**(-19.0)
    Coulomb_A = 10.0**(-10.0)
    Coulomb_factor = Coulomb_k*Coulomb_e*Coulomb_e/Coulomb_A
    Coulomb_factor = Coulomb_factor/Coulomb_e

    eV2Ha = 0.0367502
    A2Bohr = 1.88973
    eV_A_2_Ha_Bohr = eV2Ha/ A2Bohr

    eV2kcalmol = 23.0609 

    natom = len(Q)
    forces_QEq = np.zeros(shape=(natom,3))
    E = 0.0
    for i in range(natom):
        for j in range(i+1,natom):
            Rij = dist_matrix[i,j]
            X = (Rij**3.0)
            v_xij = min_dist_PBC(xyz[i,:], xyz[j,:], cell_3)
            #print (v_xij)
            dRij_dxi = v_xij / Rij
            dX_dxi = (Rij**2.0) * dRij_dxi
            Fi = Coulomb_factor*Q[i]*Q[j] * (X**(-4.0/3.0)) * dX_dxi * eV_A_2_Ha_Bohr
            E += Coulomb_factor*Q[i]*Q[j] * (X**(-1.0/3.0)) * eV2kcalmol
            forces_QEq[i,:] +=  np.array(Fi)
            forces_QEq[j,:] += -np.array(Fi)
    return forces_QEq, E

# compute distance using PBC, working only for orthorhombic cell?
def distance(x0, x1, cell):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * cell, delta - cell, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

def min_dist_PBC(x0, x1, cell):
    delta = np.array(x0 - x1)
    delta = np.where(delta >  0.5 * cell, delta - cell, delta)
    delta = np.where(delta < -0.5 * cell, delta + cell, delta)
    return delta

f  = open(file_xyz ,"rt")
f2  = open("q_all.dat" ,"w")
f3  = open("QEq.xyzfes" ,"w")

nconf = 1
while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break
    natom = int(tmp)

    tmp = f.readline().split()
    if (option_cell=="cell_3"):
        cell_3 = [float(tmp[i]) for i in range(3)]
    elif (option_cell=="cell_9"):
        cell_3 = [float(tmp[4*i]) for i in range(3)]
    elif (option_cell=="NON_ORTHO"):
        cell_3 = [float(tmp[4*i+1]) for i in range(3)]
    cell_3 = np.array(cell_3)

    atomList  = []
    atomIndex = []
    xyz = np.zeros(shape=(natom,3))
    for i in range(natom):
        tmp = f.readline().split()
        atomList.append(tmp[0])
        xyz[i,0] =  float(tmp[1])
        xyz[i,1] =  float(tmp[2])
        xyz[i,2] =  float(tmp[3])

    dist_matrix = np.array([])
    for i in range(natom):
        tmp_dist = distance(xyz, xyz[i,:], cell_3)
        dist_matrix = np.append(dist_matrix, tmp_dist)

    dist_matrix = dist_matrix.reshape((natom, natom))

    forces_QEq, E = compute_forces_QEq(q, dist_matrix, xyz, cell_3)

    f3.write("%d\n" %( natom ))
    for i in range(3):
        f3.write("%12.6f" %( cell_3[i] ))
    f3.write("%20.6f" %( E ))
    f3.write("\n")


    for i in range(natom):
        f3.write("%s" %( atomList[i] ))
        for j in range(3):
            f3.write("%12.6f" %( xyz[i,j] ))
        for j in range(3):
            f3.write("%12.6f" %( forces_QEq[i,j]  ))
        f3.write("\n")


    nconf += 1



