import numpy as np

def read_masses(filename):
    """
    Read the second column from the Masses section of a LAMMPS data file.

    Parameters
    ----------
    filename : str
        Path to the LAMMPS data file

    Returns
    -------
    np.ndarray
        Array of masses (second column)
    """
    masses = []
    in_masses = False
    data_started = False

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Enter Masses section
            if line == "Masses":
                in_masses = True
                continue

            # Skip empty line right after "Masses"
            if in_masses and not data_started and line == "":
                continue

            # Start reading data
            if in_masses and line != "":
                data_started = True
                parts = line.split()
                masses.append(float(parts[1]))
                continue

            # Stop when data section ends (empty line after data)
            if in_masses and data_started and line == "":
                break

    return np.array(masses)

def read_atoms_type_charge(filename):
    """
    Read columns 3 and 4 (type, charge) and xyz (5,6,7) from the Atoms section
    of a LAMMPS data file.

    Returns
    -------
    np.ndarray, np.ndarray: [atom_type, atom_charge]
    """
    atom_types = []
    atom_charges = []
    xyz = []
    in_atoms = False
    data_started = False

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Enter Atoms section
            if line == "Atoms":
                in_atoms = True
                continue

            # Skip empty line right after "Atoms"
            if in_atoms and not data_started and line == "":
                continue

            # Read data lines
            if in_atoms and line != "":
                data_started = True
                parts = line.split()

                # column 3: atom type
                atom_types.append(int(parts[2]))
                # column 4: charge
                atom_charges.append(float(parts[3]))
                # column 5,6,7: xyz
                xyz.append([float(parts[4]),float(parts[5]),float(parts[6])])
                continue

            # Stop at empty line after data
            if in_atoms and data_started and line == "":
                break
    return np.array(atom_types), np.array(atom_charges), np.array(xyz)

def read_bonds_type_index(filename):
    """
    Read columns 2 and [3,4] (type, index) from the Bonds section
    of a LAMMPS data file.

    Returns
    -------
    np.ndarray, np.ndarray: [bonds_type, bonds_index]
    """
    bond_types = []
    bond_index = []
    in_bonds = False
    data_started = False

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Enter Atoms section
            if line == "Bonds":
                in_bonds = True
                continue

            # Skip empty line right after "Bonds"
            if in_bonds and not data_started and line == "":
                continue

            # Read data lines
            if in_bonds and line != "":
                data_started = True
                parts = line.split()

                # column 2: bond type
                bond_types.append(int(parts[1]))
                # column 3,4: bond index
                bond_index.append([int(parts[2]),int(parts[3])])
                continue

            # Stop at empty line after data
            if in_bonds and data_started and line == "":
                break
    return np.array(bond_types), np.array(bond_index)

def read_angles_type_index(filename):
    """
    Read columns 2 and [3,4] (type, index) from the Angles section
    of a LAMMPS data file.

    Returns
    -------
    np.ndarray, np.ndarray: [angles_type, angles_index]
    """
    angle_types = []
    angle_index = []
    in_angles = False
    data_started = False

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Enter Atoms section
            if line == "Angles":
                in_angles = True
                continue

            # Skip empty line right after "Angles"
            if in_angles and not data_started and line == "":
                continue

            # Read data lines
            if in_angles and line != "":
                data_started = True
                parts = line.split()

                # column 2: angle type
                angle_types.append(int(parts[1]))
                # column 3,4,5: angle index
                angle_index.append([int(parts[2]),int(parts[3]),int(parts[4])])
                continue

            # Stop at empty line after data
            if in_angles and data_started and line == "":
                break
    return np.array(angle_types), np.array(angle_index)

def read_dihedrals_type_index(filename):
    """
    Read columns 2 and [3,4,5,6] (type, index) from the Dihedrals section
    of a LAMMPS data file.

    Returns
    -------
    np.ndarray, np.ndarray: [dihedrals_type, dihedrals_index]
    """
    dihedral_types = []
    dihedral_index = []
    in_dihedrals = False
    data_started = False

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Enter Atoms section
            if line == "Dihedrals":
                in_dihedrals = True
                continue

            # Skip empty line right after "Dihedrals"
            if in_dihedrals and not data_started and line == "":
                continue

            # Read data lines
            if in_dihedrals and line != "":
                data_started = True
                parts = line.split()

                # column 2: dihedral type
                dihedral_types.append(int(parts[1]))
                # column 3,4,5,6: dihedral index
                dihedral_index.append([int(parts[2]),int(parts[3]),int(parts[4]),int(parts[5])])
                continue

            # Stop at empty line after data
            if in_dihedrals and data_started and line == "":
                break
    return np.array(dihedral_types), np.array(dihedral_index)

def read_impropers_type_index(filename):
    """
    Read columns 2 and [3,4,5,6] (type, index) from the Impropers section
    of a LAMMPS data file.

    Returns
    -------
    np.ndarray, np.ndarray: [impropers_type, impropers_index]
    """
    improper_types = []
    improper_index = []
    in_impropers = False
    data_started = False

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Enter Atoms section
            if line == "Impropers":
                in_impropers = True
                continue

            # Skip empty line right after "Impropers"
            if in_impropers and not data_started and line == "":
                continue

            # Read data lines
            if in_impropers and line != "":
                data_started = True
                parts = line.split()

                # column 2: improper type
                improper_types.append(int(parts[1]))
                # column 3,4,5,6: improper index
                improper_index.append([int(parts[2]),int(parts[3]),int(parts[4]),int(parts[5])])
                continue

            # Stop at empty line after data
            if in_impropers and data_started and line == "":
                break
    return np.array(improper_types), np.array(improper_index)

def read_coeffs_bond(file_path):
    """
    Read lines starting with 'bond_coeff' and return:
      - bond_ids: np.array of column 2 (integers)
      - bond_params: np.array of columns 3 and 4 (floats, shape (N, 2))
    """
    bond_ids = []
    bond_params = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            # skip empty lines and pure comment lines
            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if parts[0] == "bond_coeff":
                # column 2: bond type id
                bond_ids.append(int(parts[1]))
                # columns 3 and 4: parameters
                bond_params.append([float(parts[2]), float(parts[3])])

    bond_ids = np.array(bond_ids, dtype=int)
    bond_params = np.array(bond_params, dtype=float)

    return bond_ids, bond_params

def read_coeffs_angle(file_path):
    """
    Read lines starting with 'angle_coeff' and return:
      - angle_ids: np.array of column 2 (integers)
      - angle_params: np.array of columns 3 and 4 (floats, shape (N, 2))
    """
    angle_ids = []
    angle_params = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            # skip empty lines and pure comment lines
            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if parts[0] == "angle_coeff":
                # column 2: angle type id
                angle_ids.append(int(parts[1]))
                # columns 3 and 4: parameters
                angle_params.append([float(parts[2]), float(parts[3])])

    angle_ids = np.array(angle_ids, dtype=int)
    angle_params = np.array(angle_params, dtype=float)

    return angle_ids, angle_params

def read_coeffs_dihedral(file_path):
    """
    Read lines starting with 'dihedral_coeff' and return:
      - dihedral_ids: np.array of column 2 (integers)
      - dihedral_sym: np.array of column 3 (str)
      - dihedral_params: np.array of columns 4, 5, 6 and 7 (floats, shape (N, 4))
    """
    dihedral_ids = []
    dihedral_sym = []
    dihedral_params = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            # skip empty lines and pure comment lines
            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if parts[0] == "dihedral_coeff":
                # column 2: dihedral type id
                dihedral_ids.append(int(parts[1]))
                # column 3: dihedral type symbol
                dihedral_sym.append(parts[2])
                # columns 4, 5, 6, and 7: parameters
                dihedral_params.append([float(parts[3]), float(parts[4]), float(parts[5]), float(parts[6]) ])

    dihedral_ids = np.array(dihedral_ids, dtype=int)
    dihedral_sym = np.array(dihedral_sym, dtype=str)
    dihedral_params = np.array(dihedral_params, dtype=float)

    return dihedral_ids, dihedral_sym, dihedral_params

def read_coeffs_improper(file_path):
    """
    Read lines starting with 'improper_coeff' and return:
      - improper_ids: np.array of column 2 (integers)
      - improper_params: np.array of columns 3 and 4 (floats, shape (N, 2))
    """
    improper_ids = []
    improper_params = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            # skip empty lines and pure comment lines
            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if parts[0] == "improper_coeff":
                # column 2: improper type id
                improper_ids.append(int(parts[1]))
                # columns 3 and 4: parameters
                improper_params.append([float(parts[2]), float(parts[3])])

    improper_ids = np.array(improper_ids, dtype=int)
    improper_params = np.array(improper_params, dtype=float)

    return improper_ids, improper_params

def have_same_unique_values(a, b):
    """
    Return True if arrays a and b have exactly the same unique values.
    """
    uniq_a = np.unique(a)
    uniq_b = np.unique(b)
    return np.array_equal(np.sort(uniq_a), np.sort(uniq_b))

df = {}
moles = ["water","pentanol"]
for mole in moles:
    file = f"pools/{mole}/data_{mole}.dat"
    print (file)

    masses = read_masses(file)

    atom_types, atom_charges, xyz = read_atoms_type_charge(file)
    if len(masses) > len(atom_types):
        masses = masses[atom_types-1]

    atom_type_min = np.min(atom_types)
    atom_types = atom_types-atom_type_min+1

    bond_types, bond_index = read_bonds_type_index(file)
    bond_ids, bond_params = read_coeffs_bond(file)

    if not have_same_unique_values(bond_ids,bond_types):
        print(f"Error: unique values do not match between bond_types")
        print (bond_types)
        print (bond_ids)
        sys.exit(1)

    bond_type_min = np.min(bond_types)
    bond_types = bond_types-bond_type_min+1
    bond_ids = bond_ids-bond_type_min+1

    angle_types, angle_index = read_angles_type_index(file)
    angle_ids, angle_params = read_coeffs_angle(file)

    if not have_same_unique_values(angle_ids,angle_types):
        print(f"Error: unique values do not match between angle_types")
        print (angle_types)
        print (angle_ids)
        sys.exit(1)

    angle_type_min = np.min(angle_types)
    angle_types = angle_types-angle_type_min+1
    angle_ids = angle_ids-angle_type_min+1

    dihedral_types, dihedral_index = read_dihedrals_type_index(file)
    dihedral_ids, dihedral_sym, dihedral_params = read_coeffs_dihedral(file)

    if not have_same_unique_values(dihedral_ids,dihedral_types):
        print(f"Error: unique values do not match between dihedral_types")
        print (dihedral_types)
        print (dihedral_ids)
        sys.exit(1)

    if len(dihedral_types)>0:
        dihedral_type_min = np.min(dihedral_types)
        dihedral_types = dihedral_types-dihedral_type_min+1
        dihedral_ids = dihedral_ids-dihedral_type_min+1

    improper_types, improper_index = read_impropers_type_index(file)
    improper_ids, improper_params = read_coeffs_improper(file)

    if not have_same_unique_values(improper_ids,improper_types):
        print(f"Error: unique values do not match between improper_types")
        print (improper_types)
        print (improper_ids)
        sys.exit(1)

    if len(improper_types)>0:
        improper_type_min = np.min(improper_types)
        improper_types = improper_types-improper_type_min+1
        improper_ids = improper_ids-improper_type_min+1

    df[mole] = {
        'masses':          masses,
        'atom_types':      atom_types,
        'atom_charges':    atom_charges,
        'xyz':             xyz,
        'bond_types':      bond_types,
        'bond_index':      bond_index,
        'bond_ids':        bond_ids,
        'bond_params':     bond_params,
        'angle_types':     angle_types,
        'angle_index':     angle_index,
        'angle_ids':       angle_ids,
        'angle_params':    angle_params,
        'dihedral_types':  dihedral_types,
        'dihedral_index':  dihedral_index,
        'dihedral_ids':    dihedral_ids,
        'dihedral_sym':    dihedral_sym,
        'dihedral_params': dihedral_params,
        'improper_types':  improper_types,
        'improper_index':  improper_index,
        'improper_ids':    improper_ids,
        'improper_params': improper_params
    }


print (df)
