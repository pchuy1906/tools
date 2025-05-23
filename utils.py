import numpy as np

conv_Ha_2_kcalmol = 627.503
conv_eV_2_kcalmol = 23.0609 
conv_eV_2_Ha      = 0.0367502

conv_Angstrom_2_Bohr = 1.8897259886


def calculate_cell_parameters(matrix):
    """
    Calculates cell lengths (a, b, c) and angles (alpha, beta, gamma) from a 3x3 matrix.

    Args:
        matrix (numpy.ndarray): A 3x3 matrix representing the unit cell.

    Returns:
        tuple: A tuple containing cell lengths (a, b, c) and angles (alpha, beta, gamma) in degrees.
    """
    a = np.linalg.norm(matrix[0])
    b = np.linalg.norm(matrix[1])
    c = np.linalg.norm(matrix[2])

    alpha = np.degrees(np.arccos(np.dot(matrix[1], matrix[2]) / (b * c)))
    beta = np.degrees(np.arccos(npible_float(np.dot(matrix[0], matrix[2]) / (a * c))))
    gamma = np.degrees(np.arccos(np.dot(matrix[0], matrix[1]) / (a * b)))

    return (a, b, c), (alpha, beta, gamma)

def Cell_XYZ_ABC(cellXYZ):
    a1= cellXYZ[0,:]
    a2= cellXYZ[1,:]
    a3= cellXYZ[2,:]
    a = np.sqrt(np.dot(a1, a1))
    b = np.sqrt(np.dot(a2, a2))
    c = np.sqrt(np.dot(a3, a3))
    alp = np.arccos(np.dot(a2, a3)/(b*c))*180.0/np.pi
    bet = np.arccos(np.dot(a1, a3)/(a*c))*180.0/np.pi
    gam = np.arccos(np.dot(a1, a2)/(a*b))*180.0/np.pi
    return np.array([a,b,c]), np.array([alp,bet,gam])

def create_matrix_from_lengths_angles(lengths, angles_degrees):
    """
    Constructs a 3x3 matrix from given lengths and angles (in degrees).

    Args:
        lengths: A list or numpy array of three numbers representing the lengths of the matrix's vectors.
        angles_degrees: A list or numpy array of three numbers representing the angles (in degrees) 
                        between the matrix's vectors.

    Returns:
        A 3x3 numpy array representing the constructed matrix.
    """
    angles_radians = np.radians(angles_degrees)
    alpha, beta, gamma = angles_radians
    a, b, c = lengths

    # Construct the matrix
    matrix = np.array([
        [a, b * np.cos(gamma), c * np.cos(beta)],
        [0, b * np.sin(gamma), c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)],
        [0, 0, c * np.sqrt(1 - np.cos(beta)**2 - ((np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma))**2)]
    ])

    return matrix.T


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


def read_multicharge(file_out_multicharge, natom):

    f = open(file_out_multicharge, 'rt')

    while True:

        line = f.readline()
        if line == '': break

        keywords = " q "
        if keywords in line:
            print ("*****")
            print ("read atomic charges")
            print (line)
            line = f.readline()
            q = []
            for i in range(natom):
                line = f.readline().split()
                q.append(float(line[4]))

        keywords = "dE/dx"
        if keywords in line:
            print ("*****")
            print ("read atomic forces")
            print (line)
            line = f.readline()
            fxyz = np.array([])
            for i in range(natom):
                line = f.readline().split()
                tmp = [-float(line[i]) for i in range(3,6)]
                fxyz = np.append(fxyz, np.array(tmp))
                # Unit of forces Ha/Bohr, ChIMES uses this unit
            fxyz = fxyz.reshape((natom, 3))


        keywords = "Electrostatic energy:"
        if keywords in line:
            print ("*****")
            print ("read energy")
            print (line)
            line = line.split()
            energy = float(line[2])*conv_Ha_2_kcalmol
            # Unit of energy Ha, ChIMES uses kcal/mol; need conv_Ha_2_kcalmol

        keywords = "Virial:"
        if keywords in line:
            print ("*****")
            print ("read stress")
            print (line)
            line = f.readline()
            line = f.readline()
            line = f.readline()
            stress = np.array([])
            for i in range(3):
                line = f.readline().split()
                tmp = [float(line[i]) for i in range(1,4)]
                stress = np.append(stress, np.array(tmp))
            # Unit of stress , ChIMES uses GPa
            stress = stress.reshape((3, 3))


    return q, energy, fxyz, stress


def write_xyzf(fname, natom, AtomList, xyz, cell_xyz, cell_type, fxyz, energy, stress, export_stress):
    f2 = open(fname, "w")
    f2.write("%1d\n" %( natom ))
    # cell parameters
    if (cell_type == "cell_3"):
        for i in range(3):
            f2.write("%15.9f" %( cell_xyz[i,i]))
    elif (cell_type == "cell_9"):
        for i in range(3):
            for j in range(3):
                f2.write("%15.9f" %( cell_xyz[i,j]))
    elif (cell_type == "NON_ORTHO"):
        f2.write("NON_ORTHO ")
        for i in range(3):
            for j in range(3):
                f2.write("%15.9f" %( cell_xyz[i,j]))
    else:
        print ("unknown cell_type")
    # stress
    if export_stress:
        for i in range(3):
            f2.write("%15.9f" %( stress[i,i]))
        #xy , xy, yz
        f2.write("%15.9f" %( stress[0,1]))
        f2.write("%15.9f" %( stress[0,2]))
        f2.write("%15.9f" %( stress[1,2]))
    # energy
    f2.write("%20.9f" %( energy ))
    f2.write("\n")

    for i in range(natom):
        f2.write("%s" %( AtomList[i] ))
        for j in range(3):
            f2.write("%15.9f" %( xyz[i,j]))
        for j in range(3):
            f2.write("%15.9f" %( fxyz[i,j]))
        f2.write("\n")
    
    f2.close()

