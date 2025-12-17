#!/usr/bin/env python3
"""
-----------------------------------------------------------------------------
Author: Huy Pham
Email:  pham20@llnl.gov
Description:
    Main function for bounded least-squares fitting of Lennard-Jones parameters.
-----------------------------------------------------------------------------
"""

import os
import numpy as np
import argparse
import logging

from AI4NS_utils import read_xyzf_compute_A_matrix, print_epsilon_sigma
from AI4NS_utils import gen_weights
from AI4NS_utils import lstsq_solver_linear, lstsq_solver_nonlinear


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="fit LJ parameters for xyzf file."
    )
    parser.add_argument(
        "--file_xyzf",
        default="file.xyzf",
        help="Path to the ChIMES xyzf file (default: file.xyzf)"
    )
    parser.add_argument(
        "--n_type_max",
        type=int,
        default=50,
        help="Maximum number of atom type (default: 50)"
    )
    parser.add_argument(
        "--rcut",
        type=float,
        default=10.0,
        help="Cutoff distance of atom interactions (default: 10 A)"
    )
    parser.add_argument(
        "--weight_forces",
        type=float,
        default=10.0,
        help="Weight the forces (default: 10)"
    )
    parser.add_argument(
        "--linear_fit",
        action="store_true",
        help="fit linear specified"
    )
    parser.add_argument(
        "--file_solution",
        default="x.dat",
        help="Path to the solution file (default: x.dat)"
    )
    parser.add_argument(
        "--file_weight",
        default="weight.dat",
        help="Path to the weight file (default: weight.dat)"
    )
    parser.add_argument(
        "--same_molecules",
        action="store_true",
        help="one molecule type only"
    )
    return parser.parse_args()

def main():
    """
    Main execution workflow for fitting LJ parameters
    """
    setup_logging()
    args = parse_arguments()

    file_xyzf_path = args.file_xyzf
    n_type_max     = args.n_type_max
    rcut           = args.rcut
    weight_forces  = args.weight_forces
    linear_fit     = args.linear_fit
    file_solution  = args.file_solution
    file_weight    = args.file_weight
    same_molecules = args.same_molecules

    if weight_forces > 0.0001:
        train_forces = True
    else:
        train_forces = False
    # Check if files exist
    if not os.path.isfile(file_xyzf_path):
        logging.error(f'file xyzf not found: {file_xyzf}')
        return
    logging.info('Starting to read xyzf file and generate A and b matrices')
    A_matrix, b_matrix, column_id_of, label_energy_forces, n_atom = read_xyzf_compute_A_matrix(file_xyzf_path, n_type_max, rcut, same_molecules=same_molecules, train_forces=train_forces)

    logging.info('Filtering A matrix, remove column with all zero values')
    remaining_cols = np.where(~np.all(A_matrix == 0, axis=0))[0]
    filtered_A_matrix = A_matrix[:, remaining_cols]
    inv_column_id_of = {v: k for k, v in column_id_of.items()}
    symbols_remaining_cols = [inv_column_id_of[val] for val in remaining_cols]
    print (remaining_cols)
    print (symbols_remaining_cols)
    print ("Shapes of matrix A and b:", filtered_A_matrix.shape, b_matrix.shape)
    
    if os.path.isfile(file_weight):
        weights = np.loadtxt(file_weight)
        print ("read the weights from file")
    else:
        print ("generate the energy weights based on number of atoms")
        weights = gen_weights(weight_forces, train_forces, n_atom)
    print ("Shape of the weights", len(weights))

    if linear_fit:
        logging.info('Performing constrained linear fitting')
        x = lstsq_solver_linear(filtered_A_matrix, b_matrix, weights, symbols_remaining_cols, label_energy_forces)
        direct=False
    else:
        logging.info('Performing nonlinear fitting')
        x = lstsq_solver_nonlinear(filtered_A_matrix, b_matrix, weights, symbols_remaining_cols, label_energy_forces)
        direct=True
    logging.info('Printing out pair styles for LAMMPS')
    print_epsilon_sigma(x, symbols_remaining_cols, direct=direct)

    if os.path.isfile(file_solution):
        eps_sig = np.loadtxt(file_solution)
        part_A = 4.0 * eps_sig[:,0] * eps_sig[:,1]**12
        part_B = 4.0 * eps_sig[:,0] * eps_sig[:,1]**6
        nE = len(x) - (len(part_A) + len(part_B))
        part_E = np.zeros(nE)
        x0 = np.concatenate([part_A, part_B, part_E])
        Ax = filtered_A_matrix @ x0
        np.savetxt("Ax.dat", Ax)

if __name__ == '__main__':
    main()

