from AI4NS_utils import match_mass_to_element_ase

import numpy as np
import sys
import math

def read_masses(file_path):
    """
    Read the second column from the Masses section of a LAMMPS data file.

    Parameters
    ----------
    file_path : str
        Path to the LAMMPS data file

    Returns
    -------
    np.ndarray
        Array of masses (second column)
    """
    masses = []
    in_masses = False
    data_started = False

    with open(file_path, "r") as f:
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

def read_atoms_type_charge(file_path):
    """
    Read columns 3 and 4 (type, charge) and xyz (5,6,7) from the Atoms section
    of a LAMMPS data file.

    Returns
    -------
    np.ndarray, np.ndarray: [atom_data, atom_charge]
    """
    atom_data = []
    atom_charges = []
    xyz = []
    in_atoms = False
    data_started = False

    with open(file_path, "r") as f:
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
                atom_data.append(int(parts[2]))
                # column 4: charge
                atom_charges.append(float(parts[3]))
                # column 5,6,7: xyz
                xyz.append([float(parts[4]),float(parts[5]),float(parts[6])])
                continue

            # Stop at empty line after data
            if in_atoms and data_started and line == "":
                break
    return np.array(atom_data), np.array(atom_charges), np.array(xyz)

def read_index_improper(file_path, keyword, cols):
    """
    Read columns 2 and [3,4] (type, index) from the Bonds section
    of a LAMMPS data file.

    Returns
    -------
    np.ndarray, np.ndarray: [bonds_type, bonds_index]
    """
    bond_data = []
    bond_index = []
    in_bonds = False
    data_started = False

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()

            # Enter Atoms section
            if line == keyword:
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
                bond_data.append(int(parts[1]))
                # column 3,4: bond index
                data_cols = []
                for col in cols:
                    data_cols.append(int(parts[col]))
                bond_index.append(data_cols)
                continue

            # Stop at empty line after data
            if in_bonds and data_started and line == "":
                break
    return np.array(bond_data), np.array(bond_index)

def read_coeff_improper(file_path, keyword, cols, col_sym=None):
    """
    Read lines starting with 'improper_coeff' and return:
      - improper_ids: np.array of column 2 (integers)
      - improper_params: np.array of columns 3 and 4 (floats, shape (N, 2))
    """
    improper_ids = []
    improper_params = []
    improper_sym = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            # skip empty lines and pure comment lines
            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if parts[0] == keyword:
                # column 2: improper type id
                improper_ids.append(int(parts[1]))
                # Comment
                if col_sym is not None:
                    improper_sym.append(parts[col_sym])
                # Comment
                data_cols = []
                for col in cols:
                    data_cols.append(float(parts[col]))
                improper_params.append(data_cols)

    improper_ids = np.array(improper_ids, dtype=int)
    improper_params = np.array(improper_params, dtype=float)

    if col_sym is not None:
        improper_sym = np.array(improper_sym, dtype=str)
        return improper_ids, improper_params, improper_sym
    else:
        return improper_ids, improper_params

def have_same_unique_values(a, b):
    """
    Return True if arrays a and b have exactly the same unique values.
    """
    uniq_a = np.unique(a)
    uniq_b = np.unique(b)
    return np.array_equal(np.sort(uniq_a), np.sort(uniq_b))

def gen_data_structure(moles, path_pools):
    df = {}
    for mole in moles:
        file = f"{path_pools}/data_{mole}.dat"
    
        masses = read_masses(file)
    
        atom_data, atom_charges, xyz = read_atoms_type_charge(file)
        #if len(masses) > len(atom_data):
        #    masses = masses[atom_data-1]
        masses = masses[atom_data-1]
        com = np.average(xyz, axis=0, weights=masses)
        xyz -= com

        atom_data_min = np.min(atom_data)
        atom_data = atom_data-atom_data_min+1
    
        bond_data, bond_index = read_index_improper(file,keyword='Bonds',cols=[2,3])
        bond_ids, bond_params = read_coeff_improper(file,keyword='bond_coeff',cols=[2,3])
        if not have_same_unique_values(bond_ids,bond_data):
            print(f"Error: unique values do not match between bond_data")
            print (bond_data)
            print (bond_ids)
            sys.exit(1)
        bond_data_min = np.min(bond_data)
        bond_data = bond_data-bond_data_min+1
        bond_ids = bond_ids-bond_data_min+1
    
        angle_data, angle_index = read_index_improper(file,keyword='Angles',cols=[2,3,4])
        angle_ids, angle_params = read_coeff_improper(file,keyword='angle_coeff',cols=[2,3])
        if not have_same_unique_values(angle_ids,angle_data):
            print(f"Error: unique values do not match between angle_data")
            print (angle_data)
            print (angle_ids)
            sys.exit(1)
        angle_data_min = np.min(angle_data)
        angle_data = angle_data-angle_data_min+1
        angle_ids = angle_ids-angle_data_min+1
    
        dihedral_data, dihedral_index = read_index_improper(file,keyword='Dihedrals',cols=[2,3,4,5])
        dihedral_ids, dihedral_params, dihedral_sym = read_coeff_improper(file,keyword='dihedral_coeff',cols=[3,4,5,6], col_sym=2)
        if not have_same_unique_values(dihedral_ids,dihedral_data):
            print(f"Error: unique values do not match between dihedral_data")
            print (dihedral_data)
            print (dihedral_ids)
            sys.exit(1)
        if len(dihedral_data)>0:
            dihedral_data_min = np.min(dihedral_data)
            dihedral_data = dihedral_data-dihedral_data_min+1
            dihedral_ids = dihedral_ids-dihedral_data_min+1
    
        improper_data, improper_index = read_index_improper(file,keyword='Impropers',cols=[2,3,4,5])
        improper_ids, improper_params = read_coeff_improper(file,keyword='improper_coeff',cols=[2,3])
        if not have_same_unique_values(improper_ids,improper_data):
            print(f"Error: unique values do not match between improper_data")
            print (improper_data)
            print (improper_ids)
            sys.exit(1)
        if len(improper_data)>0:
            improper_data_min = np.min(improper_data)
            improper_data = improper_data-improper_data_min+1
            improper_ids = improper_ids-improper_data_min+1
    
        df[mole] = {
            'masses':          masses,
            'atom_data':       atom_data,
            'atom_charges':    atom_charges,
            'xyz':             xyz,
            'bond_data':       bond_data,
            'bond_index':      bond_index,
            'bond_ids':        bond_ids,
            'bond_params':     bond_params,
            'angle_data':      angle_data,
            'angle_index':     angle_index,
            'angle_ids':       angle_ids,
            'angle_params':    angle_params,
            'dihedral_data':   dihedral_data,
            'dihedral_index':  dihedral_index,
            'dihedral_ids':    dihedral_ids,
            'dihedral_sym':    dihedral_sym,
            'dihedral_params': dihedral_params,
            'improper_data':   improper_data,
            'improper_index':  improper_index,
            'improper_ids':    improper_ids,
            'improper_params': improper_params
        }
    return df

def write_lmp_data_topology(f2, moles, nmoles, df, n_atom_moles, ndata, keyword1, keyword2):
    ncount_angle = 0
    ishift_type = ishift_atom = 0
    for i in range(len(moles)):
        mole = moles[i]
        nmole = nmoles[i]
        df_mole = df[mole]
        angle_data   = df_mole[keyword1]
        angle_index  = df_mole[keyword2]
        for k in range(nmole):
            for j in range(len(angle_data)):
                ncount_angle += 1
                f2.write("%d " %( ncount_angle ))
                f2.write("%d " %( angle_data[j] + ishift_type))
                for l in range(ndata):
                    f2.write("%d " %( angle_index[j,l] + ishift_atom ))
                f2.write("\n")
            ishift_atom += n_atom_moles[i]
        ishift_type += len(np.unique(angle_data))

def write_lmp_data_FFvalues(f2, moles, nmoles, df, ndata, keyword1, keyword2, keyword3=None):
    ishift_bond_type = 0
    for i in range(len(moles)):
        mole = moles[i]
        nmole = nmoles[i]
        df_mole = df[mole]
        bond_ids     = df_mole[keyword1]
        bond_params  = df_mole[keyword2]
        if keyword3 is not None:
            dihedral_sym = df_mole[keyword3]

        for j in range(len(bond_ids)):
            f2.write("%d " %( bond_ids[j] + ishift_bond_type))
            if keyword3 is not None:
                f2.write("%s " %( dihedral_sym[j]))
            for l in range(ndata):
                f2.write("%12.5f " %( bond_params[j,l] ))
            f2.write("\n")
        ishift_bond_type += len(bond_ids)

def write_lmp_data_mixture(df, moles, nmoles, maxsize):
    f2 = open("data.lammps_mixture", "w")
    n_atom = n_bond = n_angl = n_dihe = n_impr = 0
    n_atom_type = n_bond_type = n_angl_type = n_dihe_type = n_impr_type = 0

    for i in range(len(moles)):
        mole = moles[i]
        nmole = nmoles[i]
        df_mole = df[mole]

        n_atom += nmole*len(df_mole['atom_charges'])
        n_bond += nmole*len(df_mole['bond_data'])
        n_angl += nmole*len(df_mole['angle_data'])
        n_dihe += nmole*len(df_mole['dihedral_data'])
        n_impr += nmole*len(df_mole['improper_data'])

        atom_data = df_mole['atom_data']
        uniq_atom_data = np.unique(atom_data)
        n_atom_type += len(uniq_atom_data)
        n_bond_type += len(df_mole['bond_ids'])
        n_angl_type += len(df_mole['angle_ids'])
        n_dihe_type += len(df_mole['dihedral_ids'])
        n_impr_type += len(df_mole['improper_ids'])

    f2.write("%1s\n" %( "# lammps data file" ))
    f2.write("\n")
    f2.write("%1d %1s\n" %( n_atom, "atoms" ))
    f2.write("%1d %1s\n" %( n_bond, "bonds" ))
    f2.write("%1d %1s\n" %( n_angl, "angles" ))
    f2.write("%1d %1s\n" %( n_dihe, "dihedrals" ))
    f2.write("%1d %1s\n" %( n_impr, "impropers" ))
    f2.write("\n")
    f2.write("%1d %1s\n" %( n_atom_type, "atom types" ))
    f2.write("%1d %1s\n" %( n_bond_type, "bond types" ))
    f2.write("%1d %1s\n" %( n_angl_type, "angle types" ))
    f2.write("%1d %1s\n" %( n_dihe_type, "dihedral types" ))
    f2.write("%1d %1s\n" %( n_impr_type, "improper types" ))
    f2.write("\n")

    total_moles = sum(nmoles)
    nlattice = math.ceil( total_moles**(1/3) )
    all_triples = np.array(np.meshgrid(
        np.arange(nlattice), np.arange(nlattice), np.arange(nlattice)
    )).T.reshape(-1, 3)
    idx = np.random.choice(len(all_triples), size=total_moles, replace=False)
    lattice = all_triples[idx]

    latticesize = maxsize + 2.0
    f2.write("%15.6f %15.6f %1s\n" %(0.0, latticesize*nlattice, "xlo xhi" ))
    f2.write("%15.6f %15.6f %1s\n" %(0.0, latticesize*nlattice, "ylo yhi" ))
    f2.write("%15.6f %15.6f %1s\n" %(0.0, latticesize*nlattice, "zlo zhi" ))
    f2.write("%15.6f %15.6f %15.6f %1s\n" %(0.0, 0.0, 0.0, "xy xz yz" ))
    f2.write("\n")

    f2.write("%1s\n" %( "Masses" ))
    f2.write("\n")
    ishift = 0
    for i in range(len(moles)):
        mole = moles[i]
        nmole = nmoles[i]
        df_mole = df[mole]

        atom_data = df_mole['atom_data']
        masses = df_mole['masses']
        uniq_atoms, indices = np.unique(atom_data, return_index=True)
        masses_by_type = masses[indices]
        for atom_type, mass in zip(uniq_atoms, masses_by_type):
            f2.write("%1d %15.9f\n" %( atom_type+ishift, mass ))
        ishift += len(masses_by_type)
    f2.write("\n")

    f2.write("%1s\n" %( "Atoms" ))
    f2.write("\n")
    ncount_atom = ncount_mole = 0
    ishift = 0
    n_atom_moles = []
    f3 = open("config.xyz", "w")
    f3.write("%1d\n" %( n_atom ))
    f3.write("%12.6f" %( latticesize*nlattice ))
    f3.write("%12.6f" %( latticesize*nlattice ))
    f3.write("%12.6f" %( latticesize*nlattice ))
    f3.write("\n")
    for i in range(len(moles)):
        mole = moles[i]
        nmole = nmoles[i]
        df_mole = df[mole]

        atom_data    = df_mole['atom_data']
        atom_charges = df_mole['atom_charges']
        xyz          = df_mole['xyz']
        masses       = df_mole['masses']
        asym = match_mass_to_element_ase(masses, tol=0.01)
        print (asym)
        for k in range(nmole):
            R = lattice[ncount_mole,:]
            ncount_mole += 1
            for j in range(len(atom_data)):
                ncount_atom += 1
                f2.write("%d " %( ncount_atom ))
                f2.write("%d " %( ncount_mole ))
                f2.write("%d " %( atom_data[j] + ishift))
                f2.write("%12.4f " %( atom_charges[j] ))
                for l in range(3):
                    f2.write("%15.6f " %( xyz[j,l] + (R[l]+0.5) * latticesize ))
                f2.write("\n")
                f3.write("%s " %( asym[j] ))
                for l in range(3):
                    f3.write("%15.6f " %( xyz[j,l] + (R[l]+0.5) * latticesize ))
                f3.write("\n")
        ishift += len(np.unique(atom_data))
        n_atom_moles.append(len(atom_data))
    f3.close()
    f2.write("\n")

    f2.write("%1s\n" %( "Bonds" ))
    f2.write("\n")
    write_lmp_data_topology(f2, moles, nmoles, df,
        n_atom_moles,
        ndata=2, keyword1='bond_data', keyword2='bond_index',
    )
    f2.write("\n")

    f2.write("%1s\n" %( "Bond Coeffs" ))
    f2.write("\n")
    write_lmp_data_FFvalues(f2, moles, nmoles, df,
        ndata=2, keyword1='bond_ids', keyword2='bond_params',
    )
    f2.write("\n")

    f2.write("%1s\n" %( "Angles" ))
    f2.write("\n")
    write_lmp_data_topology(f2, moles, nmoles, df,
        n_atom_moles, 
        ndata=3, keyword1='angle_data', keyword2='angle_index',
    )
    f2.write("\n")

    f2.write("%1s\n" %( "Angle Coeffs" ))
    f2.write("\n")
    write_lmp_data_FFvalues(f2, moles, nmoles, df,
        ndata=2, keyword1='angle_ids', keyword2='angle_params',
    )
    f2.write("\n")

    f2.write("%1s\n" %( "Dihedrals" ))
    f2.write("\n")
    write_lmp_data_topology(f2, moles, nmoles, df,
        n_atom_moles,
        ndata=4, keyword1='dihedral_data', keyword2='dihedral_index',
    )
    f2.write("\n")

    f2.write("%1s\n" %( "Dihedral Coeffs" ))
    f2.write("\n")
    write_lmp_data_FFvalues(f2, moles, nmoles, df,
        ndata=4, keyword1='dihedral_ids', keyword2='dihedral_params',
        keyword3='dihedral_sym',
    )
    f2.write("\n")

    f2.write("%1s\n" %( "Impropers" ))
    f2.write("\n")
    write_lmp_data_topology(f2, moles, nmoles, df,
        n_atom_moles, 
        ndata=4, keyword1='improper_data', keyword2='improper_index',
    )
    f2.write("\n")

    f2.write("%1s\n" %( "Improper Coeffs" ))
    f2.write("\n")
    write_lmp_data_FFvalues(f2, moles, nmoles, df,
        ndata=2, keyword1='improper_ids', keyword2='improper_params',
    )
    f2.write("\n")


def compute_n_moles(df, moles, mole_ratio, cubic_length, density):
    Na = 6.02214076e23
    # convert from Angstrom to m
    cubic_length_in_m = cubic_length * 1e-10
    volume = cubic_length_in_m ** 3
    mole_masses = []
    l1dmax = []
    for i in range(len(moles)):
        mole = moles[i]
        df_mole = df[mole]
        masses = df_mole['masses']
        atom_data = df_mole['atom_data']
        mole_mass = np.sum(masses)
        # convert from g/mol to kg/molecule
        mole_mass *= 1e-3/Na
        mole_masses.append(mole_mass)

        xyz = df_mole['xyz']
        xmin, ymin, zmin = np.min(xyz, axis=0)
        xmax, ymax, zmax = np.max(xyz, axis=0)
        l1d = [xmax-xmin, ymax-ymin, zmax-zmin]
        l1dmax.append( max(l1d) )

    nmolA = mole_masses[0] + mole_masses[1]*(1.0-mole_ratio)/mole_ratio
    nmolA = volume * density/ nmolA
    nmolB = nmolA * (1.0 - mole_ratio)/mole_ratio
    return np.array([round(nmolA), round(nmolB)]), max(l1dmax)


moles = ["water","pentanol"]
path_pools = "./pools/"
df = gen_data_structure(moles,path_pools)

nmoles, maxsize = compute_n_moles(df, moles, mole_ratio=0.5, cubic_length=20.0, density=1000.0)
#print (nmoles)
#print (maxsize)

write_lmp_data_mixture(df, moles, nmoles, maxsize)


