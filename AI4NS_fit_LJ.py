import os
import numpy as np
import argparse
import logging

from scipy.linalg import lstsq
from scipy.linalg import svd
from scipy.optimize import lsq_linear

from AI4NS_utils import read_xyzf_compute_Amatrix


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
    Main execution workflow for fitting LF parameters
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
    Amatrix, bmatrix, = read_xyzf_compute_Amatrix(file_xyzf_path, n_type_max, rcut)

    #x, residuals, rank, s = lstsq(Amatrix, bmatrix)
    #Ax = Amatrix @ x
    #output = np.column_stack((bmatrix, Ax))
    #np.savetxt('output.dat', output, fmt='%.6f', delimiter=' ')
    #print (x)
    #print (len(x))

    remaining_cols = np.where(~np.all(Amatrix == 0, axis=0))[0]
    filtered_Amatrix = Amatrix[:, remaining_cols]
    #x, residuals, rank, s = lstsq(filtered_Amatrix, bmatrix)
    x = lsq_linear(filtered_Amatrix, bmatrix, bounds=(100.0, np.inf)).x
    Ax = filtered_Amatrix @ x
    output = np.column_stack((bmatrix, Ax))
    np.savetxt('output.dat', output, fmt='%.6f', delimiter=' ')
    print (x)
    print (len(x))

    #epsilon = 0.1
    #sigma = 5.0
    #x1 = 4.0*epsilon*sigma**12
    #x2 = 4.0*epsilon*sigma**6
    #x = np.array([x1]*18 + [x2]*18)
    #Ax = filtered_Amatrix @ x
    #output = np.column_stack((bmatrix, Ax))
    #np.savetxt('output2.dat', output, fmt='%.6f', delimiter=' ')
    #print (x)
    #print (len(x))




if __name__ == '__main__':
    main()

