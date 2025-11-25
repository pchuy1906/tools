import os
import numpy as np
import argparse
import logging
import sys
from typing import Sequence, Any
from ase.data import atomic_masses, chemical_symbols
from ase.geometry import find_mic

from numpy.linalg import lstsq
#from scipy.linalg import lstsq
from scipy.linalg import svd
from scipy.optimize import lsq_linear

from ase.units import kcal, mol, Hartree, Bohr, GPa, eV, Angstrom, bar
convert_kbar_2_GPa = bar/GPa
kcal_per_mol = kcal/mol
convert_eV_2_kcal_per_mol = eV/kcal_per_mol
convert_eV_2_Hartree = eV/Hartree
convert_Hartree_2_kcal_per_mol = Hartree/kcal_per_mol
convert_Bohr_2_Angstrom = Bohr/Angstrom
convert_Angstrom_2_Bohr = Angstrom/Bohr
convert_kcal_per_mol_2_Hartree = kcal_per_mol/Hartree

def read_poscar(filename):
    """
    Reads atomic species and their counts from a VASP POSCAR file.

    Parameters:
        filename (str): Path to the POSCAR file.

    Returns:
        tuple: (species, counts)
            species (list of str): Atomic species.
            counts (list of int): Number of atoms for each species.
    """
    with open(filename, 'rt') as file:
        # Skip the first five lines (header information)
        for _ in range(5):
            file.readline()
        # Read atomic species
        species = file.readline().split()
        ntype = len(species)
        # Read atomic counts
        counts = [int(x) for x in file.readline().split()[:ntype]]
    return species, counts

def create_atom_list(species, counts):
    """
    Generates a list of atomic species labels, repeated according to their counts.

    Parameters:
        species (list of str): List of atomic species labels.
        counts (list of int): Corresponding number of atoms for each species.

    Returns:
        list of str: Expanded list of atomic species labels.
    """
    atom_list = []
    for label, count in zip(species, counts):
        atom_list.extend([label] * count)
    return atom_list

def read_outcar(filename, natom):
    """
    Extracts cell vectors, atomic positions, forces, energy, and stress tensor from a VASP OUTCAR file.

    Parameters:
        filename (str): Path to the OUTCAR file.
        natom (int): Number of atoms in the system.

    Returns:
        tuple:
            cell_vectors (np.ndarray): Cell lattice vectors (shape: (3, 3)).
            positions (np.ndarray): Atomic positions (shape: (natom, 3)).
            forces (np.ndarray): Atomic forces (shape: (natom, 3)), converted to Hartree/Bohr.
            energy (float): Free energy in kcal/mol.
            stresses (list of float): Stress tensor components in GPa.
    """
    positions = np.zeros((natom, 3))
    forces = np.zeros((natom, 3))
    cell_vectors = np.zeros((3, 3))
    energy = None
    stresses = None

    with open(filename, 'rt') as file:
        for line in file:
            # Extract stress tensor
            if "in kB" in line:
                tokens = line.split()
                stresses = [float(tokens[i]) * convert_kbar_2_GPa for i in range(2, 8)]
            # Extract cell lattice vectors
            elif "direct lattice vectors" in line:
                for i in range(3):
                    vector_line = next(file).split()
                    cell_vectors[i] = [float(vector_line[j]) for j in range(3)]
            # Extract atomic positions and forces
            elif "POSITION" in line:
                next(file)  # Skip header line
                for i in range(natom):
                    atom_line = next(file).split()
                    positions[i] = [float(atom_line[j]) for j in range(3)]
                    forces[i] = [float(atom_line[j]) for j in range(3, 6)]
                # Convert forces to Hartree/Bohr
                forces *= convert_eV_2_Hartree / convert_Angstrom_2_Bohr
            # Extract free energy
            elif "free  energy" in line:
                energy = float(line.split()[4]) * convert_eV_2_kcal_per_mol

    return cell_vectors, positions, forces, energy, stresses

def write_chimes_xyzf(
    atom_list,
    cell_type,
    cell_vectors,
    export_stress,
    energy,
    positions,
    forces,
    stresses=None,
    filename="VASP_2_ChIMES.xyzf",
    lammps_ids=None
):
    """
    Writes atomic configuration, cell parameters, energy, forces, and optional stress to a ChIMES-format .xyzf file.

    Parameters:
        atom_list (list of str): List of atomic species labels.
        cell_type (str): Cell type ('cell_3' or 'NON_ORTHO').
        cell_vectors (np.ndarray): Cell lattice vectors (shape: (3, 3)).
        export_stress (bool): Whether to include stress tensor in the output.
        energy (float): Total energy value.
        positions (np.ndarray): Atomic positions (shape: (n_atoms, 3)).
        forces (np.ndarray): Atomic forces (shape: (n_atoms, 3)).
        stresses (list of float, optional): Stress tensor components (length 6).
        filename (str): Output filename (default: "VASP_2_ChIMES.xyzf").
    """
    natom = len(atom_list)

    if lammps_ids is not None:
        sorted_data = lammps_ids[lammps_ids[:, 1].argsort()]
        order = sorted_data[:,0]
        atype = sorted_data[:,2]
        amol  = sorted_data[:,3]


    with open(filename, "a") as file:
        # Write number of atoms
        file.write(f"{natom}\n")

        # Write cell parameters
        if cell_type == "cell_3":
            for i in range(3):
                for j in range(3):
                    file.write(f"{cell_vectors[i, j]:15.9f}")
        elif cell_type == "NON_ORTHO":
            file.write("NON_ORTHO ")
            for i in range(3):
                for j in range(3):
                    file.write(f"{cell_vectors[i, j]:15.9f}")

        # Write stresses if requested
        if export_stress:
            if stresses is None or len(stresses) != 6:
                raise ValueError("Stresses must be provided as a list of 6 values when export_stress is True.")
            for value in stresses:
                file.write(f"{value:15.9f}")

        # Write energy
        file.write(f"{energy:20.9f}\n")

        # Write atomic data
        for iatom in range(natom):

            if lammps_ids is not None:
                i = order[iatom]
            else:
                i = iatom

            file.write(f"{atom_list[i]}")
            for j in range(3):
                file.write(f"{positions[i, j]:15.9f}")
            for j in range(3):
                file.write(f"{forces[i, j]:15.9f}")

            if lammps_ids is not None:
                file.write(f" {atype[iatom]} {amol[iatom]}")
            file.write("\n")

def extract_masses(filename):
    masses = []
    f = open(filename, 'rt')
    while True:
        line = f.readline()
        if line == '': break
        if "atom types" in line:
            ntype = int(line.split()[0])
        if "Masses" in line:
            line = f.readline() # This line is empty
            for i in range(ntype):
                mass = f.readline().split()[1]
                masses.append(float(mass))
    f.close()
    return np.array(masses)


def match_mass_to_element_ase(masses, tol=0.01):
    """
    Matches each mass in the array to its element symbol using ASE data.

    Parameters:
        masses (array-like): Array of atomic masses.
        tol (float): Tolerance for matching masses.

    Returns:
        list: List of element symbols corresponding to the masses.
    """
    element_types = []
    for mass in masses:
        found = False
        for Z, std_mass in enumerate(atomic_masses):
            if abs(mass - std_mass) < tol and chemical_symbols[Z] != 'X':
                element_types.append(chemical_symbols[Z])
                found = True
                break
        if not found:
            element_types.append('Unknown')
    return np.array(element_types)


def read_cell(cell_type, line_cell_x, line_cell_y, line_cell_z):
    """
    Reads cell parameters based on the cell type and returns a NumPy array
    representing the cell in VASP/ASE format.

    Parameters:
        cell_type (str): Type of the cell ("cell_9" or "cell_3").
        line_cell_x (list): List containing x cell bounds and tilt factors.
        line_cell_y (list): List containing y cell bounds and tilt factors.
        line_cell_z (list): List containing z cell bounds and tilt factors.

    Returns:
        np.ndarray: Array of cell parameters.
    """
    if cell_type=="cell_9":
        [xlo_bound, xhi_bound, xy] = [float(x) for x in line_cell_x]
        [ylo_bound, yhi_bound, xz] = [float(x) for x in line_cell_y]
        [zlo_bound, zhi_bound, yz] = [float(x) for x in line_cell_z]
        zlo = zlo_bound
        zhi = zhi_bound
        ylo = ylo_bound - min(0.0,yz)
        yhi = yhi_bound - max(0.0,yz)
        xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
        xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
        lx = xhi-xlo
        ly = yhi-ylo
        lz = zhi-zlo
    elif cell_type=="cell_3":
        tmp = line_cell_x
        lx = float(tmp[1])-float(tmp[0])
        tmp = line_cell_y
        ly = float(tmp[1])-float(tmp[0])
        tmp = line_cell_z
        lz = float(tmp[1])-float(tmp[0])
        xy = 0.0
        xz = 0.0
        yz = 0.0
    else:
        print ("cell_type is unknown!!!", cell_type)
        sys.exit(1)
    return np.array([lx,0,0,xy,ly,0,xz,yz,lz])

def write_POSCAR(
    nconf: int,
    cell_vectors: np.ndarray,
    atom_list: Sequence[str],
    xyz: np.ndarray,
    atype: np.ndarray,
    amol: np.ndarray,
    rmin_select: float = 2.7,
) -> None:
    nmol = np.max(amol)
    molecules_xyz, molecules_atype = identify_molecules(xyz, atype, amol, nmol, fxyz=None)
    rmin = 100.0
    for i in range(nmol):
        for j in range(i+1,nmol):
            if len(molecules_atype[i]) <= len(molecules_atype[j]):
                xyz1 = molecules_xyz[i]
                atype1 = molecules_atype[i]
                xyz2 = molecules_xyz[j]
                atype2 = molecules_atype[j]
            else:
                xyz2 = molecules_xyz[i]
                atype2 = molecules_atype[i]
                xyz1 = molecules_xyz[j]
                atype1 = molecules_atype[j]
            # 
            # can we speed up here???
            for k in range(len(atype1)):
                diff = xyz2 - xyz1[k,:]
                mic_diff, _ = find_mic(diff, cell_vectors, pbc=True)
                distances = np.linalg.norm(mic_diff, axis=1)
                if rmin > np.min(distances):
                    rmin = np.min(distances)
    if rmin > rmin_select:
        poscar_filename = f"POSCAR_{nconf}"
        ntype_filename = "ntype.dat"
        order_filename = "LAMMPS_mapping.dat"
    
        # Get unique symbols and their counts
        syms, counts_syms = np.unique(atom_list, return_counts=True)
    
        new_order = []
        new_atype = []
        new_amol = []
        # Write POSCAR file
        with open(poscar_filename, "w") as f_poscar:
            f_poscar.write("COMMENT\n")
            f_poscar.write(f"{1.0:15.9f}\n")
            for vec in cell_vectors:
                f_poscar.write(f"{vec[0]:20.15f} {vec[1]:20.15f} {vec[2]:20.15f}\n")
            f_poscar.write(" " + " ".join(syms) + "\n")
            f_poscar.write(" " + " ".join(str(count) for count in counts_syms) + "\n")
            f_poscar.write("Direct\n")
    
            # Transform coordinates to fractional
            frac_xyz = np.dot(xyz, np.linalg.inv(cell_vectors))
    
            # Write atomic positions grouped by element
            natom = len(atom_list)
            for sym in syms:
                for k in range(natom):
                    if atom_list[k] == sym:
                        new_order.append(k)
                        new_atype.append(atype[k])
                        new_amol.append(amol[k])
                        f_poscar.write(
                            f"{frac_xyz[k,0]:20.15f} {frac_xyz[k,1]:20.15f} {frac_xyz[k,2]:20.15f}\n"
                        )
    
        # Write ntype.dat file
        with open(ntype_filename, "w") as f_ntype:
            for sym in syms:
                f_ntype.write(f"{sym}\n")
    
        # Write LAMMPS IDs file
        with open("LAMMPS_ID_order_type_mol.dat", "w") as file:
            for idx, (order, atype, amol) in enumerate(zip(new_order, new_atype, new_amol)):
                file.write(f"{idx} {order} {atype} {amol}\n")

def read_lammps_dump(filename, element_symbols, potential_energies=None, export_stress=None, lammps_units=None):
    """
    Reads a LAMMPS dump file and optionally associates potential energies.

    Parameters
    ----------
    lammps_dump_path : str
        Path to the LAMMPS dump file.
    element_symbols : list of str
        List of element symbols corresponding to atom types.
    potential_energies : array-like, optional
        Array of potential energy values, default is None.
    """
    with open(filename, "rt") as f:
        nconf = 0
        while True:
            line = f.readline()
            if not line:
                break  # End of file

            # Skip three lines (assuming header or irrelevant info)
            for _ in range(3):
                tmp = f.readline()

            natom_line = tmp.strip()
            try:
                natom = int(natom_line.split()[0])
            except (IndexError, ValueError) as e:
                logging.error(f"Failed to parse atom count: {natom_line} ({e})")
                print(f"Error: Failed to parse atom count: {natom_line}")
                sys.exit(1)

            # Read cell lines
            cell_info_line = f.readline()
            line_cell_x = f.readline().split()
            line_cell_y = f.readline().split()
            line_cell_z = f.readline().split()

            if "xy" in cell_info_line:
                cell_type = "cell_9"
            else:
                cell_type = "cell_3"

            cell9 = read_cell(cell_type, line_cell_x, line_cell_y, line_cell_z)

            # Read atom data header
            atom_header = f.readline().split()
            atom_header = np.array(atom_header) 
            try:
                idx = np.where(atom_header == 'x')[0][0]
                idy = np.where(atom_header == 'y')[0][0]
                idz = np.where(atom_header == 'z')[0][0]
            except IndexError:
                idx = np.where(atom_header == 'xu')[0][0]
                idy = np.where(atom_header == 'yu')[0][0]
                idz = np.where(atom_header == 'zu')[0][0]
            idtype  = np.where(atom_header=='type')[0][0]
            idmol   = np.where(atom_header=='mol')[0][0]

            if potential_energies is not None:
                idfx = np.where(atom_header == 'fx')[0][0]
                idfy = np.where(atom_header == 'fy')[0][0]
                idfz = np.where(atom_header == 'fz')[0][0]


            # Read atom data lines
            atom_data = []
            x  = np.zeros(shape=(natom))
            y  = np.zeros(shape=(natom))
            z  = np.zeros(shape=(natom))
            fx  = np.zeros(shape=(natom))
            fy  = np.zeros(shape=(natom))
            fz  = np.zeros(shape=(natom))
            atype  = np.zeros(shape=(natom), dtype=int)
            amol  = np.zeros(shape=(natom), dtype=int)

            for i in range(natom):
                atom_line = f.readline().split()
                idt = int(atom_line[0])
                x[idt-1]  = float(atom_line[idx-2])
                y[idt-1]  = float(atom_line[idy-2])
                z[idt-1]  = float(atom_line[idz-2])
                atype[idt-1]  = float(atom_line[idtype-2])
                amol[idt-1]  = float(atom_line[idmol-2])
                if potential_energies is not None:
                    fx[idt-1]  = float(atom_line[idfx-2])
                    fy[idt-1]  = float(atom_line[idfy-2])
                    fz[idt-1]  = float(atom_line[idfz-2])

            atom_list = element_symbols[atype-1]
            nconf += 1
            cell_vectors = cell9.reshape((3, 3))
            positions = np.column_stack((x, y, z))

            if potential_energies is not None:
                energy = potential_energies[nconf-1]
                forces = np.column_stack((fx, fy, fz))
                # To be included: read stress
                stresses = None
                aorder = [i for i in range(natom)]
                aorder = np.array(aorder)
                lammps_ids = np.column_stack((aorder, aorder, atype, amol))

                # Convert quantities to ChIMES units
                if lammps_units=="real":
                    forces *= convert_kcal_per_mol_2_Hartree/convert_Angstrom_2_Bohr

                write_chimes_xyzf(
                    atom_list=atom_list,
                    cell_type="NON_ORTHO",
                    cell_vectors=cell_vectors,
                    export_stress=export_stress,
                    energy=energy,
                    positions=positions,
                    forces=forces,
                    stresses=stresses,
                    lammps_ids=lammps_ids
                )

            else:
                write_POSCAR(nconf, cell_vectors, atom_list, positions, atype, amol)

def read_lammps_log(filename):
    """
    Reads the potential energy ('PotEng') values from a LAMMPS log file.
    To be included: read stress

    Parameters
    ----------
    filename : str
        Path to the LAMMPS log file.

    Returns
    -------
    np.ndarray
        Array of potential energy values extracted from the log file.
    """
    potential_energies = []
    try:
        with open(filename, "rt") as file:
            while True:
                line = file.readline()
                if not line:
                    break  # End of file
                if "Step" in line:
                    thermo_header = line.split()
                    try:
                        poteng_index = thermo_header.index('PotEng')
                    except ValueError:
                        raise ValueError("'PotEng' not found in thermo header.")
                    try:
                        thermo_values = file.readline().split()
                        potential_energies.append(float(thermo_values[poteng_index]))
                    except ValueError:
                        thermo_values = file.readline().split()
                        thermo_values = file.readline().split()
                        potential_energies.append(float(thermo_values[poteng_index]))
                if "units" in line:
                    lammps_units = line.split()[1]

    except Exception as e:
        print(f"Error reading file '{filename}': {e}")
        return np.array([])
    return lammps_units, np.array(potential_energies)

def extract_cell_xyzf(line):
    if line[0]=="NON_ORTHO":
        cell_vectors = [float(line[i]) for i in range(1,10)]
        cell_vectors = np.array(cell_vectors)
        cell_vectors = cell_vectors.reshape((3, 3))
        if len(line)==11:
            energy = float(line[10])
    return cell_vectors, energy

def mapping_ij_to_column(n_type_max):
    ncount = 0
    column_id = []
    column_str = []
    for i in range(n_type_max):
        for j in range(i+1,n_type_max):
            tmp_str = "A_"+str(i+1)+"_"+str(j+1)
            column_id.append(ncount)
            column_str.append(tmp_str)
            ncount += 1
    for i in range(n_type_max):
        for j in range(i+1,n_type_max):
            tmp_str = "B_"+str(i+1)+"_"+str(j+1)
            column_id.append(ncount)
            column_str.append(tmp_str)
            ncount += 1
    for i in range(n_type_max):
        tmp_str = "E_"+str(i+1)
        column_id.append(ncount)
        column_str.append(tmp_str)
        ncount += 1
    column_id_of = dict(zip(column_str, column_id))
    return column_id_of

def identify_molecules(xyz, atype, amol):
    molecules_xyz = []
    molecules_atype = []
    nmol = np.max(amol)
    for i in range(nmol):
        indices = np.where(amol == i+1)
        xyz0  = xyz[indices]
        atype0 = atype[indices]
        molecules_xyz.append(xyz0)
        molecules_atype.append(atype0)
    return molecules_xyz, molecules_atype

def gen_matrix_for_single_config(xyz, atype, amol, chem_formular, column_id_of, cell_vectors, rcut, n_type_max, energy, fxyz=None):



    # A_matrix has 3 block in the column
    # First block is for part with interaction A/r^12, length nA
    # Second block is for part with interaction -B/r^6, length nB
    # Third block is for atom energies, length n_type_max
    nA = n_type_max*(n_type_max-1)//2
    nB = nA
    n_cols = nA + nB + n_type_max

    # In the row, the order is first energy, then forces of corresponding molecules
    n_rows = 1 + 3*len(xyz)

    # The b_matrix is easy to generate:
    b_matrix = []
    b_matrix.append(energy)
    if fxyz is not None:
        for i in range(len(fxyz)):
            b_matrix.append(fxyz[i,0])
            b_matrix.append(fxyz[i,1])
            b_matrix.append(fxyz[i,2])
    b_matrix = np.array(b_matrix)

    A_matrix = np.zeros((n_rows, n_cols))
    A_matrix[0,-n_type_max:] = chem_formular

    # decompose into molecules for easy compute interaction energy/forces
    molecules_xyz, molecules_atype = identify_molecules(xyz, atype, amol)

    id_atom = 0
    nmol = len(molecules_xyz)
    for i in range(nmol):
        for j in range(nmol):
            # if two molecules are not the same, do the fit
            if not np.array_equal(molecules_atype[i], molecules_atype[j]):
                xyz1 = molecules_xyz[i]
                atype1 = molecules_atype[i]
                xyz2 = molecules_xyz[j]
                atype2 = molecules_atype[j]
                for k in range(len(atype1)):
                    diff = xyz1[k,:] - xyz2
                    mic_diff, _ = find_mic(diff, cell_vectors, pbc=True)
                    distances = np.linalg.norm(mic_diff, axis=1)
                    indexes = np.where(distances < rcut)[0]
                    for good_i in indexes:
                        r = distances[good_i]
                        gen_str = f"A_{atype2[good_i]}_{atype1[k]}" if atype2[good_i] < atype1[k] else f"A_{atype1[k]}_{atype2[good_i]}"
                        column_id = column_id_of[gen_str]
                        A_matrix[0,column_id] += 1.0/r**12
                        A_matrix[0,column_id+nA] += -1.0/r**6
                        A_matrix[3*id_atom+1,column_id] += 12.0/r**13 * mic_diff[good_i][0]/r
                        A_matrix[3*id_atom+2,column_id] += 12.0/r**13 * mic_diff[good_i][1]/r
                        A_matrix[3*id_atom+3,column_id] += 12.0/r**13 * mic_diff[good_i][2]/r
                        A_matrix[3*id_atom+1,column_id+nA] += -6.0/r**7 * mic_diff[good_i][0]/r
                        A_matrix[3*id_atom+2,column_id+nA] += -6.0/r**7 * mic_diff[good_i][1]/r
                        A_matrix[3*id_atom+3,column_id+nA] += -6.0/r**7 * mic_diff[good_i][2]/r
                    id_atom += 1
    return A_matrix, b_matrix

def read_xyzf_compute_A_matrix(filename, n_type_max, rcut, train_forces=False):
    try:
        A_matrix = []
        b_matrix = []
        column_id_of = mapping_ij_to_column(n_type_max)
        with open(filename, "rt") as file:
            while True:
                line = file.readline()
                if not line:
                    break  # End of file
                n_atom = int(line.split()[0])
                line = file.readline().split()
                cell_vectors, energy = extract_cell_xyzf(line)
                x, y, z, fx, fy, fz, atype, amol = ([] for _ in range(8))
                for i in range(n_atom):
                    line = file.readline().split()
                    x.append(float(line[1]))
                    y.append(float(line[2]))
                    z.append(float(line[3]))
                    fx.append(float(line[4]))
                    fy.append(float(line[5]))
                    fz.append(float(line[6]))
                    atype.append(int(line[7]))
                    amol.append(int(line[8]))
                xyz   = np.column_stack((x, y, z))
                fxyz   = np.column_stack((fx, fy, fz)) * convert_Hartree_2_kcal_per_mol/convert_Bohr_2_Angstrom
                atype = np.array(atype)
                counts = np.bincount(atype, minlength=n_type_max+1)
                chem_formular = counts[1:n_type_max+1]
                chem_formular = np.zeros(n_type_max)
                amol  = np.array(amol)

                A_sub, b_sub = gen_matrix_for_single_config(xyz, atype, amol, chem_formular, column_id_of, cell_vectors, rcut, n_type_max, energy, fxyz)

                A_matrix.append(A_sub)
                b_matrix.append(b_sub)

        A_matrix = np.vstack(A_matrix)
        b_matrix = np.concatenate(b_matrix)
    except Exception as e:
        print(f"Error reading file '{filename}': {e}")
        return np.array([]), np.array([]), np.array([])
    return A_matrix, b_matrix, column_id_of


def lstsq_solver(A_matrix, b_matrix, weights):
    #W = np.sqrt(weights)
    #A_weighted = A_matrix * W[:, np.newaxis]
    #b_weighted = b_matrix * W
    #res = lsq_linear(A_weighted, b_weighted, bounds=(0.0, np.inf))
    #res = lsq_linear(A_matrix, b_matrix, bounds=(60.0, np.inf))
    res = lsq_linear(A_matrix, b_matrix, bounds=(200, 500000))
    x = res.x

    # Find x using least squares
    #x, residuals, rank, s = lstsq(A_matrix, b_matrix)
    #x, residuals, rank, s = lstsq(A_matrix, b_matrix, rcond=None)
    print (x)
    print (len(x))
    Ax = A_matrix @ x
    rmse = np.sqrt(np.mean((Ax - b_matrix)**2))
    print ("rmse=",rmse)
    output = np.column_stack((b_matrix, Ax))
    np.savetxt('parity.dat', output, fmt='%.6f', delimiter=' ')
    return x


def print_epsilon_sigma(x, symbols_remaining_cols):
    if len(x) != len(symbols_remaining_cols):
        print ("error")
        sys.exit(1)
    else:
        # Find all CC codes for A and B
        A_codes = [s[2:] for s in symbols_remaining_cols if s.startswith("A_")]
        B_codes = [s[2:] for s in symbols_remaining_cols if s.startswith("B_")]
        # Find the intersection of CC codes present in both A and B
        common_codes = set(A_codes) & set(B_codes)
        # For each CC, find the index of "A_$CC" and "B_$CC"
        indices = []
        for cc in sorted(common_codes):
            a_label = f"A_{cc}"
            b_label = f"B_{cc}"
            a_idx = symbols_remaining_cols.index(a_label)
            b_idx = symbols_remaining_cols.index(b_label)
            An = x[a_idx] ** (1/6)
            Bn = x[b_idx]
            sigma = (An/Bn) ** (1/6)
            epsilon = Bn/(4.0*(An/Bn))
            #print (cc, epsilon, sigma)



