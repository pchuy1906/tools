import os
import numpy as np

from ase.units import kcal, mol, Hartree, Bohr, GPa, eV, Angstrom, bar
kcal_per_mol = kcal/mol
conv_Hartree_2_kcal_per_mol = Hartree/kcal_per_mol
conv_eV_2_kcal_per_mol = eV/kcal_per_mol
conv_eV_2_Hartreertree = eV/Hartree
conv_Angstrom_2_Bohr = Angstrom/Bohr
eV_per_Angstrom3 = eV/Angstrom**3
con_eV_per_Angstrom3_to_bar = eV_per_Angstrom3/bar
con_eV_per_Angstrom3_to_GPa = eV_per_Angstrom3/GPa


def read_SPARC_out(file_SPARC_out):
    """
    Reads lattice vectors from a SPARC.out file.

    Parameters
    ----------
    file_SPARC_out : str
        Path to SPARC.out file

    Returns
    -------
    np.ndarray
        3x3 numpy array with lattice vectors (in Bohr)
    """
    lattice = []
    reading = False

    with open(file_SPARC_out, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith("Lattice vectors"):
                reading = True
                continue

            if reading:
                if not line or line.startswith("Volume") or line.startswith("*"):
                    break  # stop after 3 lines or before next section
                parts = line.split()
                if len(parts) == 3:
                    lattice.append([float(x) for x in parts])

    return np.array(lattice)/conv_Angstrom_2_Bohr


import re

def read_SPARC_static(file_SPARC_static):
    atomList = []
    coords = []
    forces = []
    stress = []
    energy = None
    reading_coords = False
    reading_forces = False
    reading_stress = False

    with open(file_SPARC_static, "r") as f:
        for line in f:
            line = line.strip()

            # Detect start of coordinate block
            if line.startswith("Fractional coordinates of"):
                reading_coords = True
                element = line.split()[-1].strip(":")
                continue
            if reading_coords:
                if line == "" or line.startswith("Total free energy") or line.startswith("Fractional coordinates"):
                    reading_coords = False
                else:
                    atomList.append(element)
                    coords.append([float(x) for x in line.split()])
                    continue

            # Extract energy
            if line.startswith("Total free energy"):
                energy = float(line.split(":")[1].strip())
                continue

            # Detect start of forces
            if line.startswith("Atomic forces"):
                reading_forces = True
                continue
            if reading_forces:
                if line == "" or line.startswith("Stress"):
                    reading_forces = False
                else:
                    forces.append([float(x) for x in line.split()])
                    continue

            # Detect stress
            if line.startswith("Stress"):
                reading_stress = True
                continue
            if reading_stress:
                if line == "":
                    reading_stress = False
                else:
                    stress.append([-float(x) for x in line.split()])

    coords = np.array(coords)
    forces = np.array(forces)
    stress = np.array(stress)

    return atomList, coords, energy, forces, stress


def read_dftb_detail(file_dftb_detail, natom):
    f = open(file_dftb_detail, 'rt')
    while True:
        line = f.readline()
        if line == '': break
        keywords = "Total Forces"
        if keywords in line:
            print ("*****")
            print ("read atomic forces")
            print (line)
            fxyz = np.array([])
            for i in range(natom):
                line = f.readline().split()
                tmp = [float(line[i]) for i in range(1,4)]
                fxyz = np.append(fxyz, np.array(tmp))
                # Unit of forces Ha/Bohr, ChIMES uses this unit
        keywords = "Force related energy:"
        if keywords in line:
            print ("*****")
            print ("read energy")
            print (line)
            eV2kcalmol =  23.0609
            energy = float( line.split()[-2] ) * eV2kcalmol
        keywords = "Total stress tensor"
        if keywords in line:
            print ("*****")
            print ("read stresses")
            print (line)
            stress = np.array([])
            for i in range(3):
                line = f.readline().split()
                tmp = [float(line[i]) for i in range(3)]
                stress = np.append(stress, np.array(tmp))
            au2GPa = 1.0/(0.339893208050290 * 0.0001)
            stress = stress * au2GPa
    fxyz = fxyz.reshape((natom,3))
    stress = stress.reshape((3,3))
    f.close()
    return fxyz, energy, stress

def read_dftb_gen(file_dftb_gen):

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

    return natom, AtomList, xyz, cell_xyz

def write_input_SPARC(cell, MESH_SPACING, ECUT, DFT_METHOD):
    # Format each element to 12 decimal places
    latvec_lines = []
    for row in cell:
        line = " ".join(f"{elem:.12f}" for elem in row)
        latvec_lines.append(line)
    latvec_text = "\n".join(latvec_lines)

    text = f"""LATVEC:
{latvec_text}
LATVEC_SCALE: {conv_Angstrom_2_Bohr} {conv_Angstrom_2_Bohr} {conv_Angstrom_2_Bohr}
"""
    if MESH_SPACING:
        text += f"MESH_SPACING: {MESH_SPACING}\n"
    if ECUT:
        text += f"ECUT: {ECUT}\n"

    text += f"""FD_ORDER: 12
BC: P P P
EXCHANGE_CORRELATION: {DFT_METHOD}
ELEC_TEMP_TYPE: fermi-dirac
ELEC_TEMP: 300
TOL_SCF: 1e-6

CALC_PRES: 1
CALC_STRESS: 1
PRINT_FORCES: 1
PRINT_ATOMS: 1
"""

    file_inpt = 'SPARC.inpt'
    with open(file_inpt, 'w') as file:
        file.write(text)



def write_ion_SPARC(natom, atomList, xyz, pathPP):
    comment = (
            """
#=========================
# format of ion file
#=========================
# ATOM_TYPE:   <atom type name> 
# PSEUDO_POT:  <path/to/pseudopotential/>
# N_TYPE_ATOM: <num of atoms of this type>
# ATOMIC_MASS: <mass of atom of this type> #(optional, for MD only)
# COORD:
# <xcoord> <ycoord> <zcoord>
# ...
""")
    file_ion = 'SPARC.ion'
    with open(file_ion, 'w') as file:
        file.write(comment)
        file.write("\n")
        for i in range(len(atomList)):
            file.write("\n")
            text = (f'ATOM_TYPE: {atomList[i]}\n')
            file.write(text)
            text = (f'ATOMIC_MASS: 1.0\n')
            file.write(text)
            text = (f'N_TYPE_ATOM: {natom[i]}\n')
            file.write(text)
            files = os.listdir(pathPP)
            keywords = (f'{atomList[i]}_')
            PP_file = [os.path.join(pathPP, f) for f in files if keywords in f]
            if (len(PP_file) !=1):
                print (PP_file[0])
                print (PP_file[1])
                raise Exception("Found more than 1 PseudoPotential")
                exit()
            text = (f'PSEUDO_POT: {PP_file[0]}\n')
            file.write(text)
            text = (f'COORD:\n')
            file.write(text)
            xyz[i] *= conv_Angstrom_2_Bohr
            for row in xyz[i]:
                line = " ".join(f"{elem:.6f}" for elem in row)
                file.write(line + "\n")


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
    beta  = np.degrees(np.arccos(np.dot(matrix[0], matrix[2]) / (a * c)))
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

    tmp = f.readline().split()
    atomNumList = [int(tmp[i]) for i in range(len(tmp))]

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

    return natom, unitcell, atomList, Axyz

def read_POSCAR_SPARC(file_POSCAR):

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

    tmp = f.readline().split()
    atomNumList = [int(tmp[i]) for i in range(len(tmp))]
    index_start = []
    index_end   = []
    for i in range(len(atomNumList)):
        index_start.append(sum(atomNumList[:i]))
        index_end.append(sum(atomNumList[:i+1]))

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

    xyz_final = []
    for i in range(len(atomNumList)):
        xyz_type_i = Axyz[index_start[i]:index_end[i],:]
        xyz_final.append(xyz_type_i)

    return atomNumList, unitcell, atomNameList, xyz_final


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

def write_POSCAR(natom, cell_3_3, atomList, xyz, edit_cell=False):
    f2 = open("POSCAR", "w")
    f2.write("%1s\n" %( "COMMENT" ))
    f2.write("%15.9f\n" %( 1.0 ))

    if edit_cell:
        lengths, angles_degrees = calculate_cell_parameters(cell_3_3)
        ncell_3_3 = create_matrix_from_lengths_angles(lengths, angles_degrees)
    else:
        ncell_3_3 = cell_3_3

    f2.write("%20.15f %20.15f %20.15f\n" %( ncell_3_3[0,0], ncell_3_3[0,1], ncell_3_3[0,2] ))
    f2.write("%20.15f %20.15f %20.15f\n" %( ncell_3_3[1,0], ncell_3_3[1,1], ncell_3_3[1,2] ))
    f2.write("%20.15f %20.15f %20.15f\n" %( ncell_3_3[2,0], ncell_3_3[2,1], ncell_3_3[2,2] ))

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
            energy = float(line[2])*conv_Hartree_2_kcal_per_mol
            # Unit of energy Ha, ChIMES uses kcal/mol; need conv_Hartree_2_kcal_per_mol

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
        #xy, xz, yz
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

def read_energy_stress(file_LAMMPS_log="log.lammps"):
    f = open(file_LAMMPS_log, 'rt')
    while True:
        line = f.readline()
        if line == '': break
        keywords = "Time"
        if keywords in line:
            tmp = line.split()
            tmp = np.array(tmp)
            # search for energy
            idE = np.where(tmp=='PotEng')[0][0]
            # search for stress
            stress_keys = ['Pxx', 'Pxy', 'Pxz', 'Pxy', 'Pyy', 'Pyz', 'Pxz', 'Pyz', 'Pzz']
            stress_ids = [np.where(tmp==tmpkey)[0][0] for tmpkey in stress_keys]
            tmp_quantities = f.readline().split()
            all_quantities = np.array([float(i) for i in tmp_quantities])
            return all_quantities[idE], all_quantities[stress_ids]


def _chimes_write_xyzf(fname, atomlist, xyz, cell_xyz, fxyz, energy, stress):
    """
    Function to write fitting data to a xyz file for ChIMES LSQ.
    """
    f2 = open(fname, 'a')
    natom = len(atomlist)
    f2.write("%1d\n" % (natom))
    # cell parameters
    f2.write("NON_ORTHO ")
    for i in range(3):
        for j in range(3):
            f2.write("%9.4f" % (cell_xyz[i, j]))
    if (len(stress) > 0):
        # Voigt notation for stress tensor:
        # xx, yy, zz
        for i in range(3):
            f2.write("%12.4f" % (stress[i]))
        # stress off-diagonal xy, xz, yz
        f2.write("%12.4f" % (stress[5]))
        f2.write("%12.4f" % (stress[4]))
        f2.write("%12.4f" % (stress[3]))
    # energy
    f2.write("%20.4f" % (energy))
    f2.write("\n")
    # xyz, fxyz
    for i in range(natom):
        f2.write("%s" % (atomlist[i]))
        for j in range(3):
            f2.write("%15.9f" % (xyz[i, j]))
        for j in range(3):
            f2.write("%15.9f" % (fxyz[i, j]))
        f2.write("\n")
    f2.close()
