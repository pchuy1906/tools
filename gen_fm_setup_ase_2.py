import numpy as np
from ase import Atoms
from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii
from ase.io import read, write
from ase.geometry import  get_distances

def apair(atom_list, atom_type):
    """
    Calculates the minimum distance between any two points from two lists.

    Args:
        atom_list: A list of atom symbols: ['C', 'H', 'N', 'O'].
        atom_type: A specific atom symbol: "N".

    Returns:
        atom pairs that alphabetically sorted:['CN', 'HN', 'NN', 'NO'].
    """
    natom = len(atom_list)
    pair = []
    for i in range(natom):
        tmp_0 = [atom_type, atom_list[i]]
        tmp_0.sort()
        pair.append(tmp_0[0]+ tmp_0[1])
    return pair

def rmin_calc(list_dist, list_pair, atypes):
    """
    Calculates the minimum distance between any two points from two lists.

    Args:
        list_dist: A list of distances: ['1.1', '1.3', '1.2', '1.5'].
        list_pair: A list of distances: ['CC', 'CH', 'HH', 'CN'].
        atypes: atom types that alphabetically sorted:['C', 'H', 'N', 'O'].

    Returns:
        The minimum distance for different pairs: 
            ['CC', 'HH', 'NN', 'OO', 'CH', 'CN', 'CO', 'HN', 'HO', 'NO']
    """

    ntype = len(atypes)
    rmins = []
    # pairs with one atom type
    for i in range(ntype):
        tpair = atypes[i] + atypes[i]
        # if the pair is found, then find the minimun distance; else set minimum=100
        iloc = [j for j in range(len(list_pair)) if list_pair[j]==tpair]
        if len(iloc)==0:
            rmin = 100.0
        else:
            rmin = np.min(list_dist[iloc])
        rmins.append(rmin)
    # pairs with two atom types
    for i in range(ntype):
        for k in range(i+1,ntype):
            tpair = atypes[i] + atypes[k]
            # if the pair is found, then find the minimun distance; else set minimum=100
            iloc = [j for j in range(len(list_pair)) if list_pair[j]==tpair]
            if len(iloc)==0:
                rmin = 100.0
            else:
                rmin = np.min(list_dist[iloc])
            rmins.append(rmin)
    return rmins


# read input 
import argparse
parser = argparse.ArgumentParser(description='generate file fm_setup.in')
# Arguments supported by the code.
parser.add_argument("--file_xyz", default='file.xyz', help='file with format xyzf, xyzfe, xyzfes')
parser.add_argument('--atom_types', nargs='+')
parser.add_argument('--polynomial_orders', nargs='+', type=int)
parser.add_argument('--cutoff_distances', nargs='+', type=float)
parser.add_argument('--delta_penalty', default=0.01,type=float)


args    = parser.parse_args()
file_xyz          = args.file_xyz
atom_types        = args.atom_types
polynomial_orders = args.polynomial_orders
cutoff_distances  = args.cutoff_distances
delta_penalty     = args.delta_penalty


# print out input
print ("atom type: %s" % atom_types)
print ("Read file xyz: %s" % file_xyz)
f  = open(file_xyz ,"rt")

nconf = 0
ncondensed = 0
ntype = len(atom_types)
npair = ntype*(ntype+1)//2
rmin_1 = [100.0]* npair

Amatrix = np.array([])

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp)
    
    cell_xyz = np.zeros(shape=(3,3))
    force = np.zeros(shape=(3*natom))
    stress = np.zeros(shape=(6))
    energy = 1000

    tmp = f.readline().split()
    #print (tmp)
    if tmp[0]=="NON_ORTHO":
        is_cell_3 = False
        if len(tmp)==11:
            # format: ["NON_ORTHO", cell[0,:], cell[1,:], cell[2,:], energy]
            cell_xyz[0,:] = [float(x) for x in tmp[1:4]]
            cell_xyz[1,:] = [float(x) for x in tmp[4:7]]
            cell_xyz[2,:] = [float(x) for x in tmp[7:10]]
            energy = float(tmp[10])
            stress = 1000
            #print ("cell_xyz", cell_xyz)
        elif len(tmp)==17:
            # format: ["NON_ORTHO", cell[0,:], cell[1,:], cell[2,:]
            #           sigma_xx/yy/zz/xy/yz/zx, energy]
            cell_xyz[0,:] = [float(x) for x in tmp[1:4]]
            cell_xyz[1,:] = [float(x) for x in tmp[4:7]]
            cell_xyz[2,:] = [float(x) for x in tmp[7:10]]
            energy = float(tmp[16])
            stress = [float(x) for x in tmp[10:16]]
            ncondensed += 1
        else:
            exit()
    else:
        is_cell_3 = True
        if len(tmp)==4:
            # format: [a,b,c, energy]
            for i in range(3):
                cell_xyz[i,i] = float(tmp[i])
            energy = float(tmp[3])
            stress = 1000
            cell_3 = [float(tmp[i]) for i in range(3)]

        elif len(tmp)==10:
            # format: [a, b, c, sigma_xx/yy/zz/xy/yz/zx, energy]
            for i in range(3):
                cell_xyz[i,i] = float(tmp[i])
            energy = float(tmp[9])
            stress = [float(x) for x in tmp[3:9]]
            ncondensed += 1
            cell_3 = [float(tmp[i]) for i in range(3)]
        else:
            exit()

    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(natom):
        tmp = f.readline().split()
        atomList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
        force[3*k+0] = float(tmp[4])
        force[3*k+1] = float(tmp[5])
        force[3*k+2] = float(tmp[6])

    atoms = Atoms(symbols=atomList, positions=xyz, cell=cell_xyz, pbc=True)
    print (atoms)
    cell = atoms.get_cell()

    nconf += 1

    pair = []
    dist = np.array([])
    for i in range(natom):
        tmp_arr = get_distances(atoms.positions[i], atoms.positions,pbc=True, cell=cell)
        tmp_dist = tmp_arr[1][0]
        tmp_dist[tmp_dist<0.01] = 100.0
        tmp_pair = apair(atomList, atomList[i])
        pair.extend(tmp_pair)
        dist = np.append(dist, tmp_dist)
    #print (len(dist), len(pair))
    rmin_2 = rmin_calc(dist, pair, atom_types)
    Amatrix = np.append(Amatrix, rmin_2)
    #print 
    #print ("minimum distances in the frame %d" %nconf)
    #print (rmin_2)
    rmin_1 = np.minimum(rmin_1, rmin_2)

if (nconf*npair != len(Amatrix)):
    print ("ERROR")
    exit()
Amatrix = Amatrix.reshape((nconf, npair))
np.savetxt('rmin.dat', Amatrix, fmt = '%.6f')

f.close
print 
print ("Total number of configuration is %s" % nconf)
print ("number of condensed phase is %s" % ncondensed)
print
print ("minimum distances in training set:")
print (rmin_1)
np.savetxt('rmin_all.dat', rmin_1, fmt = '%.6f')


def write_input(file_xyz, atom_types, rmin, polynomial_orders, cutoff_distances):
    f2 = open('_fm_setup.in', "w")
    f2.write("\n")
    f2.write("####### CONTROL VARIABLES #######\n")
    f2.write("\n")
    f2.write("# TRJFILE #\n")
    f2.write("%s\n" %file_xyz)
    f2.write("# WRAPTRJ #\n")
    f2.write("true\n")
    f2.write("# NFRAMES #\n")
    f2.write("%d\n" %nconf)
    f2.write("# NLAYERS #\n")
    f2.write("1\n")
    f2.write("# FITSTRS #\n")
    if ncondensed>0:
        f2.write("FIRSTALL %d\n" %ncondensed)
    else:
        f2.write("false\n")
    f2.write("# FITENER #\n")
    f2.write("true\n")
    f2.write("# FITCOUL #\n")
    f2.write("false\n")
    f2.write("# FITPOVR #\n")
    f2.write("false\n")
    f2.write("# PAIRTYP #\n")
    f2.write("CHEBYSHEV %d %d %d\n" %(polynomial_orders[0], polynomial_orders[1], polynomial_orders[2]) )
    f2.write("# CHBTYPE #\n")
    f2.write("MORSE\n")
    f2.write("# SPLITFI #\n")
    f2.write("true\n")
    f2.write("# SKIP_FRAMES #\n")
    f2.write("1\n")
    f2.write("\n")
    f2.write("####### TOPOLOGY VARIABLES #######\n")
    f2.write("\n")
    f2.write("# NATMTYP # \n")
    f2.write("%d\n" %len(atom_types))
    f2.write("\n")
    f2.write("# TYPEIDX #     # ATM_TYP #     # ATMCHRG #     # ATMMASS #\n")
    for i in range(len(atom_types)):
        atomic_number = atomic_numbers[atom_types[i]]
        atomic_mass = atomic_masses[atomic_number]
        f2.write("%4d %4s %4d %8.3f\n" %(i+1, atom_types[i], 0, atomic_mass ))
    f2.write("\n")
    f2.write("# PAIRIDX # ")
    f2.write("# ATM_TY1 # ")    
    f2.write("# ATM_TY1 # ")
    f2.write("# S_MINIM # ")
    f2.write("# S_MAXIM # ")
    f2.write("# S_DELTA # ")
    f2.write("# MORSE_LAMBDA # ")
    f2.write("# USEOVRP # ")    
    f2.write("# NIJBINS # ")    
    f2.write("# NIKBINS # ")    
    f2.write("# NJKBINS #\n")

    ncount = 0
    for i in range(len(atom_types)):
        ncount += 1
        tpair = atom_types[i]+atom_types[i]
        atomic_number_1 = atomic_numbers[atom_types[i]]
        atomic_number_2 = atomic_numbers[atom_types[i]]
        atomic_radius_1 = covalent_radii[atomic_number_1]
        atomic_radius_2 = covalent_radii[atomic_number_2]
        bond_length = atomic_radius_1 + atomic_radius_2
        f2.write("%4s " %(ncount))
        f2.write("%4s " %(atom_types[i]))
        f2.write("%4s " %(atom_types[i]))
        f2.write("%7.3f " %(rmin[ncount-1]))
        f2.write("%7.3f " %(cutoff_distances[0]))
        f2.write("%7.3f " %(0.1))
        f2.write("%7.3f " %(bond_length))
        f2.write("%s " %('false'))
        f2.write("%d %d %d \n" %(0,0,0))
    for i in range(len(atom_types)):
        for j in range(i+1,len(atom_types)):
            ncount += 1
            tpair = atom_types[i]+atom_types[j]
            atomic_number_1 = atomic_numbers[atom_types[i]]
            atomic_number_2 = atomic_numbers[atom_types[j]]
            atomic_radius_1 = covalent_radii[atomic_number_1]
            atomic_radius_2 = covalent_radii[atomic_number_2]
            bond_length = atomic_radius_1 + atomic_radius_2
            f2.write("%4s " %(ncount))
            f2.write("%4s " %(atom_types[i]))
            f2.write("%4s " %(atom_types[j]))
            f2.write("%7.3f " %(rmin[ncount-1]))
            f2.write("%7.3f " %(cutoff_distances[0]))
            f2.write("%7.3f " %(0.1))
            f2.write("%7.3f " %(bond_length))
            f2.write("%s " %('false'))
            f2.write("%d %d %d \n" %(0,0,0))

    f2.write("\n")
    f2.write("SPECIAL 3B S_MAXIM: ALL %7.3f\n" %cutoff_distances[1])
    f2.write("SPECIAL 4B S_MAXIM: ALL %7.3f\n" %cutoff_distances[2])
    f2.write("\n")
    f2.write("# FCUTTYP #\n")
    f2.write("TERSOFF 0.95\n")
    f2.write("\n")
    f2.write("# ENDFILE #\n")

write_input(file_xyz, atom_types, rmin_1-delta_penalty, polynomial_orders, cutoff_distances)


