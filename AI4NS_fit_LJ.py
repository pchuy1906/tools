import os
import numpy as np
import argparse
import logging

from AI4NS_utils import read_xyzf


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
    return parser.parse_args()

def main():
    """
    Main execution workflow for fitting LF parameters
    """
    setup_logging()
    args = parse_arguments()

    file_xyzf_path = args.file_xyzf
    n_type_max = args.n_type_max

    # Check if files exist
    if not os.path.isfile(file_xyzf_path):
        logging.error(f'file xyzf not found: {file_xyzf}')
        return
    A,b = read_xyzf(file_xyzf_path, n_type_max)

if __name__ == '__main__':
    main()

