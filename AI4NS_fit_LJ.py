import os
import numpy as np
import argparse
import logging

#from scipy.linalg import lstsq
#from scipy.linalg import svd
#from scipy.optimize import lsq_linear

from AI4NS_utils import read_xyzf_compute_A_matrix, lstsq_solver, print_epsilon_sigma


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

    # Check if files exist
    if not os.path.isfile(file_xyzf_path):
        logging.error(f'file xyzf not found: {file_xyzf}')
        return
    A_matrix, b_matrix, column_id_of = read_xyzf_compute_A_matrix(file_xyzf_path, n_type_max, rcut, train_forces=True)

    # Filter the A_matrix, remove zero columns
    remaining_cols = np.where(~np.all(A_matrix == 0, axis=0))[0]
    filtered_A_matrix = A_matrix[:, remaining_cols]
    inv_column_id_of = {v: k for k, v in column_id_of.items()}
    symbols_remaining_cols = [inv_column_id_of[val] for val in remaining_cols]
    
    weights = np.ones(len(b_matrix))
    x = lstsq_solver(filtered_A_matrix, b_matrix, weights)
    print_epsilon_sigma(x, symbols_remaining_cols)

    epsilon = 0.1
    sigma = 5.0
    x1 = 4.0*epsilon*sigma**12
    x2 = 4.0*epsilon*sigma**6
    x = np.array([x1]*18 + [x2]*18)
    Ax = filtered_A_matrix @ x
    output = np.column_stack((b_matrix, Ax))
    np.savetxt('output2.dat', output, fmt='%.6f', delimiter=' ')
    rmse = np.sqrt(np.mean((Ax - b_matrix)**2))
    print (x)
    print (len(x))
    print (rmse)



if __name__ == '__main__':
    main()

