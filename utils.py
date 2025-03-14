import numpy as np

def read_gen(file_dftb_gen):

    f = open(file_dftb_gen, 'rt')

    line = f.readline().split()
    natom = int(line[0])
    type_xyz = line[1]
    print ("number of atom = %d" %natom)

    atypes = f.readline().split()
    xyz = np.array([])
    cell_xyz = np.array([])
    AtomList = []

    for i in range(natom):

        line = f.readline().split()

        AtomList.append( atypes[ int(line[1])-1 ] )

        tmp = [float(line[j]) for j in range(2,5)]
        xyz = np.append(xyz, np.array(tmp))

    line = f.readline()

    line = f.readline().split()
    tmp = [float(line[j]) for j in range(3)]
    cell_xyz = np.append(cell_xyz, np.array(tmp))

    line = f.readline().split()
    tmp = [float(line[j]) for j in range(3)]
    cell_xyz = np.append(cell_xyz, np.array(tmp))

    line = f.readline().split()
    tmp = [float(line[j]) for j in range(3)]
    cell_xyz = np.append(cell_xyz, np.array(tmp))

    f.close()

    cell_xyz = cell_xyz.reshape((3,3))
    xyz = xyz.reshape((natom,3))
    if (type_xyz=="F"):
        xyz = np.dot(xyz,cell_xyz)

    return natom, cell_xyz, AtomList, xyz

def read_xyz(file_xyz, cell_type):
    f  = open(file_xyz ,"r")
    
    natom = f.readline()
    natom = int(natom)
    print ("the number of atom:", natom)
    
    cell_3_3 = np.zeros(shape=(3,3))
    tmp = f.readline().split()
    
    if cell_type=="cell_3":
        for k in range(3):
            cell_3_3[k,k] = float(tmp[k])
    elif cell_type=="cell_9":
        cell_3_3[0,0] = float(tmp[0])
        cell_3_3[0,1] = float(tmp[1])
        cell_3_3[0,2] = float(tmp[2])
    
        cell_3_3[1,0] = float(tmp[3])
        cell_3_3[1,1] = float(tmp[4])
        cell_3_3[1,2] = float(tmp[5])
    
        cell_3_3[2,0] = float(tmp[6])
        cell_3_3[2,1] = float(tmp[7])
        cell_3_3[2,2] = float(tmp[8])
    elif cell_type=="NON_ORTHO":
        cell_3_3[0,0] = float(tmp[1])
        cell_3_3[0,1] = float(tmp[2])
        cell_3_3[0,2] = float(tmp[3])
    
        cell_3_3[1,0] = float(tmp[4])
        cell_3_3[1,1] = float(tmp[5])
        cell_3_3[1,2] = float(tmp[6])
    
        cell_3_3[2,0] = float(tmp[7])
        cell_3_3[2,1] = float(tmp[8])
        cell_3_3[2,2] = float(tmp[9])
    else:
        print ("unknown cell_type", cell_type)
        exit()
    
    atomList = []
    xyz = np.zeros(shape=(natom,3))
    for k in range(natom):
        tmp = f.readline()
        tmp = tmp.split()
        atomList.append(tmp[0])
        xyz[k,0] =  float(tmp[1])
        xyz[k,1] =  float(tmp[2])
        xyz[k,2] =  float(tmp[3])
    f.close
    return natom, cell_3_3, atomList, xyz


def read_POSCAR(file_POSCAR):

    f  = open(file_POSCAR ,"r")
    tmp = f.readline()
    tmp = f.readline()

    unitcell = np.zeros(shape=(3,3))
    tmp = f.readline().split()
    unitcell[0,:] = tmp
    tmp = f.readline().split()
    unitcell[1,:] = tmp
    tmp = f.readline().split()
    unitcell[2,:] = tmp

    tmp = f.readline().split()
    atomNameList = [tmp[i] for i in range(len(tmp))]
    #print ("atomNameList=",atomNameList)

    tmp = f.readline().split()
    atomNumList = [int(tmp[i]) for i in range(len(tmp))]
    #print ("atomNumList=",atomNumList)

    atomList = []
    for i in range(len(atomNumList)):
        for j in range(atomNumList[i]):
            atomList.append(atomNameList[i])
    natom = sum(atomNumList)

    xyztype = f.readline().split()[0]

    xyz = np.zeros(shape=(natom,3))
    for k in range(natom):
        tmp = f.readline().split()
        xyz[k,:] =  tmp[:3]
    f.close()

    if xyztype=="direct" or xyztype=="Direct":
        Axyz = np.dot(xyz, unitcell)
    elif xyztype=="Cartesian":
        Axyz = xyz
    else:
        print ("unknown option xyztype", xyztype)
        exit()

    #print ("natom=", natom)
    #print ("unitcell=", unitcell)
    #print ("atomList=", atomList)
    #print ("xyz=", Axyz)
    return natom, unitcell, atomList, Axyz


def write_gen(fname, natom, unitcell, atomList, xyz):
    f2 = open(fname, "w")
    f2.write("%-d %4s\n" %(natom, "S"))
    asym_list = list(set(atomList))
    asym_list = np.unique(atomList)
    for i in range(len(asym_list)):
        f2.write("%s " %(asym_list[i]))
    f2.write("\n")

    
    for i in range(natom):
        id_sym = np.where(asym_list == atomList[i])[0][0]
        f2.write("%-5d %5d %15.9f %15.9f %15.9f\n" %(i+1, id_sym+1, xyz[i,0], xyz[i,1], xyz[i,2] ))
    f2.write("%15.9f%15.9f%15.9f\n" % (0.0,0.0,0.0))
    for ixyz in range(3):
        for icellxyz in unitcell[ixyz,:]:
            f2.write("%15.9f" % icellxyz)
        f2.write("\n")
    f2.close()

def write_xyz(fname, natom, cell_type, cell, atomList, xyz):
    f2 = open( fname, "w")
    f2.write("%1d\n" %( natom ))
    if cell_type == "cell_9":
        for i in range(3):
            for j in range(3):
                f2.write("%15.9f " %( cell[i][j] ))
        f2.write("\n")
    elif cell_type == "NON_ORTHO":
        f2.write("NON_ORTHO ")
        for i in range(3):
            for j in range(3):
                f2.write("%15.9f " %( cell[i][j] ))
        f2.write("\n")
    elif cell_type == "cell_3":
        for i in range(3):
            f2.write("%15.9f " %( cell[i][i] ))
        f2.write("\n")
    else:
        print ("unknown cell_type", cell_type)
        exit()

    for i in range(natom):
        f2.write("%s " %( atomList[i] ))
        for j in range(3):
            f2.write("%15.9f " %( xyz[i][j] ))
        f2.write("\n")

def write_POSCAR(natom, cell_3_3, atomList, xyz):
    f2 = open("POSCAR", "w")
    f2.write("%1s\n" %( "COMMENT" ))
    f2.write("%15.9f\n" %( 1.0 ))
    f2.write("%20.15f %20.15f %20.15f\n" %( cell_3_3[0,0], cell_3_3[0,1], cell_3_3[0,2] ))
    f2.write("%20.15f %20.15f %20.15f\n" %( cell_3_3[1,0], cell_3_3[1,1], cell_3_3[1,2] ))
    f2.write("%20.15f %20.15f %20.15f\n" %( cell_3_3[2,0], cell_3_3[2,1], cell_3_3[2,2] ))

    syms, counts_syms = np.unique(atomList, return_counts=True)

    print (syms)
    print (counts_syms)

    f3 = open("ntype.dat", "w")
    for item in syms:
        f2.write("%s %s" %(" ", item))
        f3.write("%s\n" %( item))

    f2.write("\n")
    for item in counts_syms:
        f2.write("%s %s" %(" ", item))
    f2.write("\n")
    f2.write("%1s\n" %( "Direct" ))


    xyz = np.dot(xyz, np.linalg.inv(cell_3_3))

    for item in syms:
        for k in range(0,natom):
            if (item == atomList[k]):
                f2.write("%20.15f %20.15f %20.15f\n" %( xyz[k,0], xyz[k,1], xyz[k,2] ))

    f2.close()

