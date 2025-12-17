import numpy as np

import argparse
parser = argparse.ArgumentParser(description='generate zero LJ parameters')
# Arguments supported by the code.
parser.add_argument('--types_mol_1', nargs='+', type=int)
parser.add_argument('--types_mol_2', nargs='+', type=int)

args        = parser.parse_args()
types_mol_1 = args.types_mol_1
types_mol_2 = args.types_mol_2

if np.array_equal(types_mol_1,types_mol_2):
    nlen = len(types_mol_1)
    for i in range(nlen):
        for j in range(i,nlen):
            epsilon = 0.0
            sigma   = 0.0
            print(f"pair_coeff {types_mol_1[i]:4d} {types_mol_1[j]:4d} lj/cut {epsilon:12.5f} {sigma:12.5f}")


