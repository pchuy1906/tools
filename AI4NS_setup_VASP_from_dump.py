#!/usr/bin/env python3
"""
Script to set up VASP calculations from LAMMPS dump and data files.
"""

import os
import sys
import argparse
import logging
import numpy as np
from typing import Sequence, Any
from ase.data import atomic_masses, chemical_symbols

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Setup VASP calculations from LAMMPS dump and data files.'
    )
    parser.add_argument(
        '--file_LAMMPS_dump',
        default='dump.lammps',
        help='Path to the LAMMPS dump file'
    )
    parser.add_argument(
        '--file_LAMMPS_data',
        default='data.lammps',
        help='Path to the LAMMPS data file'
    )
    parser.add_argument(
        '--file_LAMMPS_log',
        default='non_existence_file',
        help='Path to the LAMMPS log file'
    )
    return parser.parse_args()

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
    cell9: np.ndarray,
    atomList: Sequence[str],
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    atype: np.ndarray,
    amol: np.ndarray,
) -> None:
    poscar_filename = f"POSCAR_{nconf}"
    ntype_filename = "ntype.dat"
    order_filename = "LAMMPS_mapping.dat"

    # Reshape cell9 to (3, 3)
    if cell9.size != 9:
        logging.error("cell9 must have exactly 9 elements.")
        return
    cell_matrix = cell9.reshape((3, 3))

    # Get unique symbols and their counts
    syms, counts_syms = np.unique(atomList, return_counts=True)

    new_order = []
    new_atype = []
    new_amol = []
    # Write POSCAR file
    with open(poscar_filename, "w") as f_poscar:
        f_poscar.write("COMMENT\n")
        f_poscar.write(f"{1.0:15.9f}\n")
        for vec in cell_matrix:
            f_poscar.write(f"{vec[0]:20.15f} {vec[1]:20.15f} {vec[2]:20.15f}\n")
        f_poscar.write(" " + " ".join(syms) + "\n")
        f_poscar.write(" " + " ".join(str(count) for count in counts_syms) + "\n")
        f_poscar.write("Direct\n")

        # Transform coordinates to fractional
        xyz = np.column_stack((x, y, z))
        frac_xyz = np.dot(xyz, np.linalg.inv(cell_matrix))

        # Write atomic positions grouped by element
        natom = len(x)
        for sym in syms:
            for k in range(natom):
                if atomList[k] == sym:
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
    #with open("LAMMPS_ID_order_type_mol.dat", "w") as file:
    #    for order, atype, amol in zip(new_order, new_atype, new_amol):
    #        file.write(f"{order} {atype} {amol}\n")
    with open("LAMMPS_ID_order_type_mol.dat", "w") as file:
        for idx, (order, atype, amol) in enumerate(zip(new_order, new_atype, new_amol)):
            file.write(f"{idx} {order} {atype} {amol}\n")

def read_lammps_dump(filename, element_symbols, potential_energies=None):
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

            #logging.info(f"Config {nconf}: Atom count = {natom}")

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
            #logging.info(f"Cell parameters: {cell9}")

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

            atomList = element_symbols[atype-1]
            nconf += 1
            if potential_energies is not None:

            else:
                write_POSCAR(nconf, cell9, atomList, x, y, z, atype, amol)

def read_lammps_log(filename):
    """
    Reads the potential energy ('PotEng') values from a LAMMPS log file.

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
                    thermo_values = file.readline().split()
                    potential_energies.append(float(thermo_values[poteng_index]))
    except Exception as e:
        print(f"Error reading file '{filename}': {e}")
        return np.array([])
    return np.array(potential_energies)

def main():
    setup_logging()
    args = parse_arguments()

    lammps_dump_path = args.file_LAMMPS_dump
    lammps_data_path = args.file_LAMMPS_data
    lammps_log_path = args.file_LAMMPS_log

    logging.info('Starting collecting XYZ from LAMMPS dump format.')

    # Check if files exist
    if not os.path.isfile(lammps_data_path):
        logging.error(f'LAMMPS data file not found: {lammps_data_path}')
        return
    if not os.path.isfile(lammps_dump_path):
        logging.error(f'LAMMPS dump file not found: {lammps_dump_path}')
        return

    # Your conversion logic here
    logging.info(f'Reading masses from LAMMPS data file: {lammps_data_path}')
    masses = extract_masses(lammps_data_path)
    print(', '.join([f'{m:.3f}' for m in masses]))
    logging.info(f'Matching these masses with element symbols')
    element_symbols = match_mass_to_element_ase(masses)
    print (element_symbols)

    if not os.path.isfile(lammps_log_path):
        logging.info(f'Reading LAMMPS dump and writing POSCAR files')
        read_lammps_dump(lammps_dump_path, element_symbols)
    else:
        pe = read_lammps_log(lammps_log_path)
        read_lammps_dump(lammps_dump_path, element_symbols, potential_energies=pe)



if __name__ == '__main__':
    main()
