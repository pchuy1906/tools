#!/usr/bin/env python3
"""
-----------------------------------------------------------------------------
Author: Huy Pham
Email:  pham20@llnl.gov
Description:
    Main function for: 
        + setting up VASP calculations from LAMMPS dump file
        + export to xyzf format from LAMMPS dump/log files
-----------------------------------------------------------------------------
"""

import os
import sys
import argparse
import logging
import numpy as np
from typing import Sequence, Any
from ase.data import atomic_masses, chemical_symbols

from AI4NS_utils import extract_masses, match_mass_to_element_ase
from AI4NS_utils import read_lammps_log, read_lammps_dump

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
    parser.add_argument(
        "--export_stress",
        action="store_true",
        help="Export stress data if specified"
    )
    return parser.parse_args()

def main():
    setup_logging()
    args = parse_arguments()

    lammps_dump_path = args.file_LAMMPS_dump
    lammps_data_path = args.file_LAMMPS_data
    lammps_log_path = args.file_LAMMPS_log
    export_stress = args.export_stress

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
        lammps_units, pe = read_lammps_log(lammps_log_path)
        read_lammps_dump(lammps_dump_path, element_symbols, potential_energies=pe, export_stress=export_stress, lammps_units=lammps_units)

if __name__ == '__main__':
    main()
