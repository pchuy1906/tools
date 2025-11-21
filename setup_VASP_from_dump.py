#!/usr/bin/env python3
"""
Script to set up VASP calculations from LAMMPS dump and data files.
"""

import os
import sys
import argparse
import logging
import numpy as np

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

from ase.data import atomic_masses, chemical_symbols

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


def read_dump(filename, element_symbols):
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
            write_POSCAR(nconf, cell9, atomList, x, y, z, )

def main():
    setup_logging()
    args = parse_arguments()

    lammps_dump_path = args.file_LAMMPS_dump
    lammps_data_path = args.file_LAMMPS_data

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
    #print (masses)
    #print(np.array2string(masses, separator=', '))
    print(', '.join([f'{m:.3f}' for m in masses]))
    logging.info(f'Matching these masses with element symbols')
    element_symbols = match_mass_to_element_ase(masses)
    print (element_symbols)

    logging.info(f'Reading LAMMPS dump file: {lammps_dump_path}')
    read_dump(lammps_dump_path, element_symbols)


if __name__ == '__main__':
    main()
