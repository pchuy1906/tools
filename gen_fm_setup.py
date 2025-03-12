import numpy as np

# Dictionary of all elements matched with their atomic masses.
elements_dict_mass = {'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'NE' : 20.180, 'NA' : 22.990, 'MG' : 24.305,\
                 'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'CL' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'CA' : 40.078,\
                 'SC' : 44.956, 'TI' : 47.867, 'V' : 50.942, 'CR' : 51.996,\
                 'MN' : 54.938, 'Fe' : 55.845, 'CO' : 58.933, 'NI' : 58.693,\
                 'Cu' : 63.546, 'ZN' : 65.38, 'GA' : 69.723, 'GE' : 72.631,\
                 'AS' : 74.922, 'SE' : 78.971, 'BR' : 79.904, 'KR' : 84.798,\
                 'RB' : 84.468, 'SR' : 87.62, 'Y' : 88.906, 'ZR' : 91.224,\
                 'NB' : 92.906, 'MO' : 95.95, 'TC' : 98.907, 'RU' : 101.07,\
                 'RH' : 102.906, 'PD' : 106.42, 'AG' : 107.868, 'CD' : 112.414,\
                 'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,\
                 'I' : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,\
                 'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,\
                 'PM' : 144.913, 'SM' : 150.36, 'EU' : 151.964, 'GD' : 157.25,\
                 'TB' : 158.925, 'DY': 162.500, 'HO' : 164.930, 'ER' : 167.259,\
                 'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,\
                 'TA' : 180.948, 'W' : 183.84, 'RE' : 186.207, 'OS' : 190.23,\
                 'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,\
                 'TL' : 204.383, 'PB' : 207.2, 'Bi' : 208.980, 'PO' : 208.982,\
                 'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,\
                 'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U' : 238.029,\
                 'NP' : 237, 'PU' : 244, 'AM' : 243, 'CM' : 247, 'BK' : 247,\
                 'CT' : 251, 'ES' : 252, 'FM' : 257, 'MD' : 258, 'NO' : 259,\
                 'LR' : 262, 'RF' : 261, 'DB' : 262, 'SG' : 266, 'BH' : 264,\
                 'HS' : 269, 'MT' : 268, 'DS' : 271, 'RG' : 272, 'CN' : 285,\
                 'NH' : 284, 'FL' : 289, 'MC' : 288, 'LV' : 292, 'TS' : 294,\
                 'OG' : 294}

pair_dict_bond = {\
'CuCu' : 1.0,\
'CuC' : 1.0,\
'CuH' : 1.0,\
'CuN' : 1.0,\
'CuO' : 1.0,\
'CC' : 1.54,\
'CH' : 1.09,\
'CN' : 1.48,\
'CO' : 1.16,\
'HH' : 0.74,\
'HN' : 1.01,\
'HO' : 0.97,\
'NN' : 1.35,\
'NO' : 1.27,\
'OO' : 1.21, \
'AlAl' : 2.4,\
'AlC' : 2.0,\
'AlH' : 1.6,\
'AlN' : 1.9,\
'AlO' : 1.7,\
'OSi' : 1.63, \
'HSi' : 1.48, \
'SiSi' : 2.36, \
'BiBi' : 3.00, \
'CFe' : 3.00, \
'FeH' : 3.00, \
'FeFe' : 2.74 \
}

def distance(r1, r2, cell):
    """
    Compute the distance between two points in a triclinic box with PBC.
    Args:
        r1 (numpy.ndarray): Position vector of the first point.
        r2 (numpy.ndarray): Position vector of the second point.
        cell (numpy.ndarray): 3x3 matrix representing the unit cell vectors.
    Returns:
        float: Distance between the two points, accounting for PBC.
    """
    delta_r = r1 - r2
    delta_r_frac = np.linalg.solve(cell.T, delta_r)
    delta_r_frac_pbc = delta_r_frac - np.round(delta_r_frac)
    delta_r_pbc = np.matmul(cell.T , delta_r_frac_pbc)
    dist = np.linalg.norm(delta_r_pbc)
    #vector = r1 - r2
    #vector_frac = np.dot(vector, np.linalg.inv(unit_cell))  # Convert to fractional coordinates
    #vector_frac = vector_frac - np.round(vector_frac) # apply PBC
    #vector = np.dot(vector_frac, unit_cell) # convert back to cartesian coordinates
    #distance = np.linalg.norm(vector)
    return dist


def distance_cell3(x0, x1, cell):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * cell, delta - cell, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

# merge 2 atom symbols to pair, sorted alphabetically
def apair(atomList, atom):
    natom = len(atomList)
    pair = []
    for i in range(natom):
        tmp_0 = [atom, atomList[i]]
        tmp_0.sort()
        pair.append(tmp_0[0]+ tmp_0[1])
    return pair

# 
def rmin_calc(dist, pair, atom_types):
    ntype = len(atom_types)
    _rmin = []
    # same atom type
    for i in range(ntype):
        tpair = atom_types[i] + atom_types[i]
        #print (tpair)
        iloc = [j for j in range(len(pair)) if pair[j]==tpair]
        if len(iloc)==0:
            rmin = 100.0
        else:
            rmin = np.min(dist[iloc])
        _rmin.append(rmin)
    # 2 different atom types
    for i in range(ntype):
        for k in range(i+1,ntype):
            tpair = atom_types[i] + atom_types[k]
            #print (tpair)
            iloc = [j for j in range(len(pair)) if pair[j]==tpair]
            if len(iloc)==0:
                rmin = 100.0
            else:
                rmin = np.min(dist[iloc])
            _rmin.append(rmin)
    return _rmin


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
            # format: ["NON_ORTHO", cell[0,:], cell[1,:], cell[2,:], sigma_xx/yy/zz/xy/yz/zx, energy]
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

    nconf += 1

    pair = []
    dist = np.array([])
    for i in range(natom):
        if is_cell_3:
            cell_3 = np.array(cell_3)
            #print (cell_3)
            tmp_dist = distance_cell3(xyz, xyz[i,:], cell_3)
            tmp_dist[tmp_dist<0.01] = 100.0
        else:
            # much slower calculations
            tmp_dist = []
            for j in range(natom):
                if i==j:
                    tmp_dist.append(100.0)
                else:
                    tmp_dist.append(distance(xyz[j,:], xyz[i,:], cell_xyz))
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
    #f2.write("%-d %4s\n" %(natom, "S"))
    #f2.write("%-s \n" %(asym_list))
    f2.write("\n")
    f2.write("####### CONTROL VARIABLES #######\n")
    f2.write("\n")
    f2.write("# TRJFILE # ! The .xyzf file containing the trajectory. Like a typical xyz, but comment line has box dimes, and each line includes x,y, and z force\n")
    f2.write("%s\n" %file_xyz)
    f2.write("# WRAPTRJ # ! Does the trajectory file need wrapping? (i.e. post-run PBC)\n")
    f2.write("true\n")
    f2.write("# NFRAMES # ! How many frames are in the .xyzf file?\n")
    f2.write("%d\n" %nconf)
    f2.write("# NLAYERS # ! x,y, and z supercells.. small unit cell should have >= 1\n")
    f2.write("1\n")
    f2.write("# FITSTRS #\n")
    if ncondensed>0:
        f2.write("FIRSTALL %d\n" %ncondensed)
    else:
        f2.write("false\n")
    f2.write("# FITENER # ! Fit ENERGY\n")
    f2.write("true\n")
    f2.write("# FITCOUL # ! Fit charges? If false, use user-specified fixed charges, and subtract them from the forces -- NOTE: FUNCTIONALITY CURRENTLY ONLY SUPPORTED FOR TRUE, AND FALSE WITH CHARGES = 0\n")
    f2.write("false\n")
    f2.write("# FITPOVR #  ! Use ReaxFF linear overcoordination parameters? If this is false and parameters are provided below, those parameters will be subtracted from forces\n")
    f2.write("false\n")
    f2.write("# PAIRTYP # ! Short-range interaction type. See manual for accepted types. Case sensitive\n")
    f2.write("CHEBYSHEV %d %d %d\n" %(polynomial_orders[0], polynomial_orders[1], polynomial_orders[2]) )
    f2.write("# CHBTYPE # ! Are we transforming distance in terms of inverse distance(INVRSE_R), a morse-type function? (MORSE)? .. \"DEFAULT\" for no transformation.\n")
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
        f2.write("%4d %4s %4d %8.3f\n" %(i+1, atom_types[i], 0, elements_dict_mass.get(atom_types[i]) ))
    f2.write("\n")
    f2.write("# PAIRIDX #     # ATM_TY1 #     # ATM_TY1 #     # S_MINIM #     # S_MAXIM #     # S_DELTA #     # MORSE_LAMBDA #        # USEOVRP #     # NIJBINS #     # NIKBINS #     # NJKBINS #\n")

    ncount = 0
    for i in range(len(atom_types)):
        ncount += 1
        tpair = atom_types[i]+atom_types[i]
        print (tpair, ncount, atom_types[i], atom_types[i], rmin[ncount-1], cutoff_distances[0], 0.1, pair_dict_bond.get(tpair))
        f2.write("%4d %4s %4s %7.3f %7.3f %7.3f %7.3f %s %d %d %d\n" %(ncount, atom_types[i], atom_types[i], rmin[ncount-1], cutoff_distances[0], 0.1, pair_dict_bond.get(tpair), 'false', 0,0,0))
    for i in range(len(atom_types)):
        for j in range(i+1,len(atom_types)):
            ncount += 1
            tpair = atom_types[i]+atom_types[j]
            f2.write("%4d %4s %4s %7.3f %7.3f %7.3f %7.3f %s %d %d %d\n" %(ncount, atom_types[i], atom_types[j], rmin[ncount-1], cutoff_distances[0], 0.1, pair_dict_bond.get(tpair), 'false', 0,0,0))
    f2.write("\n")
    f2.write("SPECIAL 3B S_MAXIM: ALL %7.3f\n" %cutoff_distances[1])
    f2.write("SPECIAL 4B S_MAXIM: ALL %7.3f\n" %cutoff_distances[2])
    f2.write("\n")
    f2.write("# FCUTTYP #\n")
    f2.write("TERSOFF 0.95\n")
    f2.write("\n")
    f2.write("# ENDFILE #\n")

write_input(file_xyz, atom_types, rmin_1-delta_penalty, polynomial_orders, cutoff_distances)


