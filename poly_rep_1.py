import numpy as np

# compute distance using PBC, working only for orthorhombic cell?
def distance(x0, x1, cell):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * cell, delta - cell, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

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
npair = ntype*(ntype+1)/2
rmin_1 = [100.0]* npair

while True:
    tmp  = f.readline()
    line = tmp.strip()
    if line == '': break

    natom = int(tmp)
    
    cell = np.zeros(shape=(3,3))
    force = np.zeros(shape=(3*natom))
    stress = np.zeros(shape=(6))
    energy = 1000

    tmp = f.readline().split()
    try:
        # format: [a, b, c, sigma_xx/yy/zz/xy/yz/zx, energy]
        cell_abc = [float(x) for x in tmp[0:3]]
        stress = [float(x) for x in tmp[3:9]]
        energy = float(tmp[9])
        ncondensed += 1
    except:
        try:
            # format: ["NON_ORTHO", cell[0,:], cell[1,:], cell[2,:], sigma_xx/yy/zz/xy/yz/zx, energy]
            cell[0,:] = [float(x) for x in tmp[1:4]]
            cell[1,:] = [float(x) for x in tmp[4:7]]
            cell[2,:] = [float(x) for x in tmp[7:10]]
            stress = [float(x) for x in tmp[10:16]]
            energy = float(tmp[16])
            ncondensed += 1
            cell_abc = [ cell[0,0], cell[1,1], cell[2,2] ]
        except:
            # format: [a,b,c, energy]
            cell_abc = [float(x) for x in tmp[0:3]]
            stress = 1000
            energy = float(tmp[3])

    cell_abc = np.array(cell_abc)
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
        tmp_dist = distance(xyz, xyz[i,:], cell_abc)
        tmp_dist[tmp_dist<0.01] = 100.0
        tmp_pair = apair(atomList, atomList[i])
        pair.extend(tmp_pair)
        dist = np.append(dist, tmp_dist)
    #print (len(dist), len(pair))
    rmin_2 = rmin_calc(dist, pair, atom_types)
    print 
    print ("minimum distances in the frame %d" %nconf)
    print (rmin_2)
    rmin_1 = np.minimum(rmin_1, rmin_2)

f.close
print 
print ("Total number of configuration is %s" % nconf)
print ("number of condensed phase is %s" % ncondensed)
print
print ("minimum distances in training set:")
print (rmin_1)


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


