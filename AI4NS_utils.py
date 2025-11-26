import os
import numpy as np
import argparse
import logging
import sys

from typing import Sequence, Any, Tuple, Optional, List, Dict

from ase.data import atomic_masses, chemical_symbols
from ase.geometry import find_mic

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
    """
    Writes a POSCAR file for VASP from molecular configuration data.

    Parameters:
        nconf (int): Configuration number, used for file naming.
        cell_vectors (np.ndarray): 3x3 array of cell vectors.
        atom_list (Sequence[str]): List of atom symbols (e.g., ['H', 'O']).
        xyz (np.ndarray): Nx3 array of atomic coordinates.
        atype (np.ndarray): N array of atom types (integer indices).
        amol (np.ndarray): N array of molecule indices.
        rmin_select (float, optional): Minimum allowed distance between molecules to write POSCAR.
    
    Returns:
        None. Writes POSCAR and related files to disk.
    """
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

def extract_cell_xyzf(line: Sequence[str]) -> Tuple[np.ndarray, Optional[float]]:
    """
    Extracts cell vectors and energy from a line of data.

    Parameters:
        line (Sequence[str]): Input data line, expected to contain cell vector information.
            Format: ["NON_ORTHO", v1, v2, ..., v9, (optional: energy)]

    Returns:
        Tuple[np.ndarray, Optional[float]]: A tuple containing:
            - cell_vectors (np.ndarray): 3x3 array of cell vectors.
            - energy (float or None): Energy value if present, otherwise None.
    """
    cell_vectors = None
    energy = None

    if line[0] == "NON_ORTHO":
        try:
            cell_vectors = np.array([float(line[i]) for i in range(1, 10)]).reshape((3, 3))
        except (ValueError, IndexError) as e:
            raise ValueError(f"Error parsing cell vectors: {e}")

        if len(line) == 11:
            try:
                energy = float(line[10])
            except ValueError:
                raise ValueError(f"Error parsing energy value: {line[10]}")

    return cell_vectors, energy

def mapping_ij_to_column(n_type_max: int) -> dict:
    """
    Generates a mapping from interaction types to column indices.

    For each pair (i, j) where i < j and for each type i:
        - "A_i_j" and "B_i_j" are generated for pairs.
        - "E_i" is generated for each type.

    Parameters:
        n_type_max (int): The maximum number of types.

    Returns:
        dict: Dictionary mapping interaction type strings to column indices.
    """
    column_map = {}
    column_counter = 0

    # Generate "A_i_j" and "B_i_j" for each unique pair (i, j)
    for prefix in ["A", "B"]:
        for i in range(n_type_max):
            for j in range(i + 1, n_type_max):
                key = f"{prefix}_{i + 1}_{j + 1}"
                column_map[key] = column_counter
                column_counter += 1

    # Generate "E_i" for each type
    for i in range(n_type_max):
        key = f"E_{i + 1}"
        column_map[key] = column_counter
        column_counter += 1

    return column_map

def identify_molecules(
    xyz: np.ndarray,
    atype: np.ndarray,
    amol: np.ndarray
) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """
    Identifies and groups atomic coordinates and types by molecule.

    Parameters:
        xyz (np.ndarray): Array of atomic coordinates, shape (N, 3).
        atype (np.ndarray): Array of atomic types, shape (N,).
        amol (np.ndarray): Array of molecule indices, shape (N,).

    Returns:
        Tuple[List[np.ndarray], List[np.ndarray]]:
            - molecules_xyz: List of arrays, each containing the coordinates of atoms in a molecule.
            - molecules_atype: List of arrays, each containing the types of atoms in a molecule.
    """
    molecules_xyz = []
    molecules_atype = []
    nmol = int(np.max(amol))

    for i in range(nmol):
        indices = np.where(amol == i + 1)[0]
        xyz_molecule = xyz[indices]
        atype_molecule = atype[indices]
        molecules_xyz.append(xyz_molecule)
        molecules_atype.append(atype_molecule)

    return molecules_xyz, molecules_atype

def gen_matrix_for_single_config(
    xyz: np.ndarray,
    atype: np.ndarray,
    amol: np.ndarray,
    chem_formular: np.ndarray,
    column_id_of: Dict[str, int],
    cell_vectors: np.ndarray,
    rcut: float,
    n_type_max: int,
    energy: float,
    fxyz: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generates the A and b matrices for a single molecular configuration.

    The A matrix encodes interaction terms and atomic energies; the b matrix encodes target energy and optional forces.

    Parameters:
        xyz (np.ndarray): Atomic coordinates, shape (N_atoms, 3).
        atype (np.ndarray): Atomic types, shape (N_atoms,).
        amol (np.ndarray): Molecule indices for each atom, shape (N_atoms,).
        chem_formular (np.ndarray): Chemical formula, shape (n_type_max,).
        column_id_of (Dict[str, int]): Mapping from interaction string to column index.
        cell_vectors (np.ndarray): Cell vectors, shape (3, 3).
        rcut (float): Cutoff radius for interactions.
        n_type_max (int): Number of atom types.
        energy (float): Target energy value.
        fxyz (Optional[np.ndarray]): Atomic forces, shape (N_atoms, 3), optional.

    Returns:
        Tuple[np.ndarray, np.ndarray]: (A_matrix, b_matrix)
            - A_matrix: Matrix of interaction terms and atomic energies.
            - b_matrix: Target energy and forces.
    """
    nA = n_type_max * (n_type_max - 1) // 2
    nB = nA
    n_cols = nA + nB + n_type_max

    n_atoms = len(xyz)
    n_rows = 1 + 3 * n_atoms if fxyz is not None else 1

    # Initialize b_matrix: first element is energy, followed by force components if provided
    b_matrix = [energy]
    if fxyz is not None:
        b_matrix.extend(fxyz.flatten())
    b_matrix = np.array(b_matrix)

    # Initialize A_matrix
    A_matrix = np.zeros((n_rows, n_cols))
    # Set chemical formula in the last n_type_max columns of the first row
    A_matrix[0, -n_type_max:] = chem_formular

    # Group atoms by molecule
    molecules_xyz, molecules_atype = identify_molecules(xyz, atype, amol)

    atom_counter = 0
    nmol = len(molecules_xyz)
    for i in range(nmol):
        for j in range(nmol):
            if not np.array_equal(molecules_atype[i], molecules_atype[j]):
                xyz1 = molecules_xyz[i]
                atype1 = molecules_atype[i]
                xyz2 = molecules_xyz[j]
                atype2 = molecules_atype[j]
                for k, atom_type_k in enumerate(atype1):
                    diff = xyz1[k, :] - xyz2
                    mic_diff, _ = find_mic(diff, cell_vectors, pbc=True)
                    distances = np.linalg.norm(mic_diff, axis=1)
                    valid_indices = np.where(distances < rcut)[0]
                    for idx in valid_indices:
                        r = distances[idx]
                        atom_type_j = atype2[idx]
                        # Ensure consistent ordering for interaction string
                        if atom_type_j < atom_type_k:
                            interaction_str = f"A_{atom_type_j}_{atom_type_k}"
                        else:
                            interaction_str = f"A_{atom_type_k}_{atom_type_j}"
                        col_id = column_id_of[interaction_str]
                        # Energy terms
                        A_matrix[0, col_id] += 0.5 / r ** 12
                        A_matrix[0, col_id + nA] += -0.5 / r ** 6
                        # Force terms
                        if fxyz is not None:
                            force_prefactors = [
                                (12.0 / r ** 13, col_id),
                                (-6.0 / r ** 7, col_id + nA)
                            ]
                            for prefactor, column in force_prefactors:
                                for dim in range(3):  # x, y, z
                                    row_idx = 3 * atom_counter + 1 + dim
                                    A_matrix[row_idx, column] += prefactor * mic_diff[idx][dim] / r
                    atom_counter += 1

    return A_matrix, b_matrix


from typing import Tuple, Dict
import numpy as np

def read_xyzf_compute_A_matrix(
    filename: str,
    n_type_max: int,
    rcut: float,
    train_forces: bool = False
) -> Tuple[np.ndarray, np.ndarray, Dict[str, int]]:
    """
    Reads an XYZF file, computes A and b matrices for all configurations,
    and returns the stacked matrices along with the column mapping.

    Parameters:
        filename (str): Path to the input XYZF file.
        n_type_max (int): Maximum number of atom types.
        rcut (float): Cutoff radius for interactions.
        train_forces (bool): Whether to include atomic forces in training data.

    Returns:
        Tuple[np.ndarray, np.ndarray, Dict[str, int]]:
            - A_matrix: Stacked matrix of interaction terms and atomic energies.
            - b_matrix: Stacked target energy and forces.
            - column_id_of: Mapping from interaction string to column index.
    """
    try:
        print(
            "Energies and atomic forces are included in the training data"
            if train_forces else
            "Only energies are included in the training data"
        )

        A_matrix_list = []
        b_matrix_list = []
        column_id_of = mapping_ij_to_column(n_type_max)
        nconf = 0

        with open(filename, "rt") as file:
            while True:
                line = file.readline()
                if not line:
                    break  # End of file

                n_atom = int(line.split()[0])
                header_line = file.readline().split()
                cell_vectors, energy = extract_cell_xyzf(header_line)

                x, y, z, fx, fy, fz, atype, amol = ([] for _ in range(8))
                for _ in range(n_atom):
                    atom_line = file.readline().split()
                    x.append(float(atom_line[1]))
                    y.append(float(atom_line[2]))
                    z.append(float(atom_line[3]))
                    fx.append(float(atom_line[4]))
                    fy.append(float(atom_line[5]))
                    fz.append(float(atom_line[6]))
                    atype.append(int(atom_line[7]))
                    amol.append(int(atom_line[8]))

                xyz = np.column_stack((x, y, z))
                fxyz = (
                    np.column_stack((fx, fy, fz)) *
                    convert_Hartree_2_kcal_per_mol / convert_Bohr_2_Angstrom
                )
                atype_arr = np.array(atype)
                counts = np.bincount(atype_arr, minlength=n_type_max + 1)
                chem_formular = counts[1:n_type_max + 1]
                amol_arr = np.array(amol)

                if not train_forces:
                    fxyz = None

                A_sub, b_sub = gen_matrix_for_single_config(
                    xyz, atype_arr, amol_arr, chem_formular,
                    column_id_of, cell_vectors, rcut, n_type_max, energy, fxyz
                )

                A_matrix_list.append(A_sub)
                b_matrix_list.append(b_sub)
                nconf += 1

        A_matrix = np.vstack(A_matrix_list)
        b_matrix = np.concatenate(b_matrix_list)
        print(f"Number of configurations: {nconf}")
        print(f"Shapes of matrices A and b: {A_matrix.shape}, {b_matrix.shape}")

    except Exception as e:
        print(f"Error reading file '{filename}': {e}")
        return np.array([]), np.array([]), {}

    return A_matrix, b_matrix, column_id_of

def lstsq_solver(A_matrix, b_matrix, weights, symbols_remaining_cols):
    print (symbols_remaining_cols)
    count_A = sum(s.startswith("A_") for s in symbols_remaining_cols)
    count_B = sum(s.startswith("B_") for s in symbols_remaining_cols)
    count_E = sum(s.startswith("E_") for s in symbols_remaining_cols)
    arr_A = np.full(count_A, 0.001 * 2**12)
    arr_B = np.full(count_B, 0.001 * 2**6)
    arr_E = np.full(count_E, -500)
    lower_bounds = np.concatenate([arr_A, arr_B, arr_E])
    arr_A = np.full(count_A, 0.4 * 4**12)
    arr_B = np.full(count_B, 0.4 * 4**6)
    arr_E = np.full(count_E, 500)
    upper_bounds = np.concatenate([arr_A, arr_B, arr_E])
    W = np.sqrt(weights)
    A_weighted = A_matrix * W[:, np.newaxis]
    b_weighted = b_matrix * W
    print (lower_bounds)
    print (upper_bounds)
    res = lsq_linear(A_weighted, b_weighted, bounds=(lower_bounds, upper_bounds))
    x = res.x
    print (x[-count_E:])
    Ax = A_matrix @ x
    rmse = np.sqrt(np.mean((Ax - b_matrix)**2))
    print ("rmse=",rmse)
    output = np.column_stack((b_matrix, Ax))
    np.savetxt('parity.dat', output, fmt='%.6f', delimiter=' ')
    return x

def print_epsilon_sigma(x: List[float], symbols_remaining_cols: List[str]) -> None:
    """
    Prints Lennard-Jones epsilon and sigma parameters for atom pairs,
    based on provided coefficients and symbol labels.

    Parameters:
        x (List[float]): Coefficient values, must match symbols_remaining_cols in length.
        symbols_remaining_cols (List[str]): List of symbol strings (e.g., "A_1_2", "B_1_2").

    Output:
        Prints formatted pair_coeff lines for each atom pair with both A and B coefficients.
    """
    if len(x) != len(symbols_remaining_cols):
        print("Error: Length of coefficients does not match number of symbols.")
        sys.exit(1)

    # Extract CC codes for A and B interactions
    A_codes = [symbol[2:] for symbol in symbols_remaining_cols if symbol.startswith("A_")]
    B_codes = [symbol[2:] for symbol in symbols_remaining_cols if symbol.startswith("B_")]

    # Find common CC codes present in both A and B
    common_codes = set(A_codes) & set(B_codes)

    for cc in sorted(common_codes):
        a_label = f"A_{cc}"
        b_label = f"B_{cc}"
        try:
            a_idx = symbols_remaining_cols.index(a_label)
            b_idx = symbols_remaining_cols.index(b_label)
        except ValueError:
            print(f"Error: Could not find indices for {a_label} or {b_label}.")
            continue

        An = x[a_idx]
        Bn = x[b_idx]

        # Calculate sigma and epsilon according to Lennard-Jones parameters
        try:
            sigma = (An / Bn) ** (1 / 6)
            epsilon = Bn / (4.0 * (An / Bn))
        except ZeroDivisionError:
            print(f"Error: Division by zero for pair {cc}.")
            continue

        try:
            style_list = [int(i) for i in cc.split('_')]
        except ValueError:
            print(f"Error: Invalid CC code format '{cc}'.")
            continue

        print(f"pair_coeff {style_list[0]:4d} {style_list[1]:4d} lj/cut {epsilon:12.5f} {sigma:12.5f}")









