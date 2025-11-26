#!/usr/bin/env python3
"""
-----------------------------------------------------------------------------
Author: Huy Pham
Email:  pham20@llnl.gov
Description:
    Convert VASP to xyzf format, keeping track of LAMMPS atom types.
-----------------------------------------------------------------------------
"""

import os
import numpy as np
import argparse
import logging


from ase.units import kcal, mol, Hartree, Bohr, GPa, eV, Angstrom, bar
convert_kbar_2_GPa = bar/GPa
kcal_per_mol = kcal/mol
convert_eV_2_kcal_per_mol = eV/kcal_per_mol
convert_eV_2_Hartree = eV/Hartree
convert_Angstrom_2_Bohr = Angstrom/Bohr

from AI4NS_utils import read_poscar, read_outcar
from AI4NS_utils import create_atom_list, write_chimes_xyzf

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Convert VASP output files for ChIMES input."
    )
    parser.add_argument(
        "--vasp_outcar",
        default="OUTCAR",
        help="Path to the VASP OUTCAR file (default: OUTCAR)"
    )
    parser.add_argument(
        "--vasp_poscar",
        default="POSCAR",
        help="Path to the VASP POSCAR file (default: POSCAR)"
    )
    parser.add_argument(
        "--lammps_ids",
        default="LAMMPS_ID.dat",
        help="Path to the LAMMPS ID file (default: LAMMPS_ID.dat)"
    )
    parser.add_argument(
        "--cell_type",
        default="cell_3",
        choices=["cell_3", "NON_ORTHO"],
        help="Cell type: 'cell_3' or 'NON_ORTHO' (default: cell_3)"
    )
    parser.add_argument(
        "--export_stress",
        action="store_true",
        help="Export stress data if specified"
    )
    return parser.parse_args()

def main():
    """
    Main execution workflow for converting VASP output files to ChIMES input format.
    """
    setup_logging()
    args = parse_arguments()

    vasp_poscar_path = args.vasp_poscar
    vasp_outcar_path = args.vasp_outcar
    cell_type = args.cell_type
    export_stress = args.export_stress
    lammps_ids_path = args.lammps_ids

    # Check if files exist
    if not os.path.isfile(vasp_poscar_path):
        logging.error(f'VASP POSCAR file not found: {vasp_poscar_path}')
        return
    if not os.path.isfile(vasp_outcar_path):
        logging.error(f'VASP OUTCAR file not found: {vasp_outcar_path}')
        return

    if not os.path.isfile(lammps_ids_path):
        lammps_ids = None
    else:
        lammps_ids = np.loadtxt(lammps_ids_path, dtype=int)

    try:
        # Read atomic species and counts from POSCAR
        species, counts = read_poscar(vasp_poscar_path)
        atom_list = create_atom_list(species, counts)
        logging.info(f"Reading POSCAR")

        # Read geometry, forces, energy, and stresses from OUTCAR
        n_atoms = len(atom_list)
        cell_vectors, positions, forces, energy, stresses = read_outcar(vasp_outcar_path, n_atoms)
        logging.info(f"Reading OUTCAR")

        # Write data to ChIMES xyzf format
        write_chimes_xyzf(
            atom_list=atom_list,
            cell_type=cell_type,
            cell_vectors=cell_vectors,
            export_stress=export_stress,
            energy=energy,
            positions=positions,
            forces=forces,
            stresses=stresses,
            lammps_ids=lammps_ids
        )
        logging.info("ChIMES xyzf file successfully written.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == '__main__':
    main()

