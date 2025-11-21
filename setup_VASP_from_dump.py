#!/usr/bin/env python3
"""
Script to set up VASP calculations from LAMMPS dump and data files.
"""

import os
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
            print (line)
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
    return element_types

def main():
    setup_logging()
    args = parse_arguments()

    lammps_dump_path = args.file_LAMMPS_dump
    lammps_data_path = args.file_LAMMPS_data

    logging.info('Starting collecting XYZ from LAMMPS dump format.')
    logging.info(f'Reading LAMMPS data file: {lammps_data_path}')
    logging.info(f'Reading LAMMPS dump file: {lammps_dump_path}')

    # Check if files exist
    if not os.path.isfile(lammps_data_path):
        logging.error(f'LAMMPS data file not found: {lammps_data_path}')
        return
    if not os.path.isfile(lammps_dump_path):
        logging.error(f'LAMMPS dump file not found: {lammps_dump_path}')
        return

    # Your conversion logic here
    logging.info(f'Reading masses from LAMMPS dump file: {lammps_dump_path}')
    masses = extract_masses(lammps_data_path)
    print (masses)
    element_symbols = match_mass_to_element_ase(masses)
    print (element_symbols)

if __name__ == '__main__':
    main()
