import numpy as np

# compute distance using PBC, working only for orthorhombic cell?
def distance(x0, x1, cell):
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


args    = parser.parse_args()
file_xyz          = args.file_xyz
atom_types        = args.atom_types


# print out input
print ("atom type: %s" % atom_types)
print ("Read file xyz: %s" % file_xyz)
f  = open(file_xyz ,"rt")

nconf = 0
ncondensed = 0
ntype = len(atom_types)
npair = ntype*(ntype+1)/2
npair = int(npair)
rmin_1 = [100.0]* npair

Amatrix = np.array([])

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
        ncondensed += 1
    except:
        try:
            # format: ["NON_ORTHO", cell[0,:], cell[1,:], cell[2,:], sigma_xx/yy/zz/xy/yz/zx, energy]
            cell[0,:] = [float(x) for x in tmp[1:4]]
            cell[1,:] = [float(x) for x in tmp[4:7]]
            cell[2,:] = [float(x) for x in tmp[7:10]]
            ncondensed += 1
            cell_abc = [ cell[0,0], cell[1,1], cell[2,2] ]
        except:
            # format: [a,b,c, energy]
            cell_abc = [float(x) for x in tmp[0:3]]

    cell_abc = np.array(cell_abc)
    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(natom):
        tmp = f.readline().split()
        atomList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])

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


