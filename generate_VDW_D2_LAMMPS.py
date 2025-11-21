data = {
    "H": {"C6": 0.14, "R0": 1.001},
    "C": {"C6": 1.75, "R0": 1.452},
    "N": {"C6": 1.23, "R0": 1.397},
    "O": {"C6": 0.70, "R0": 1.342}
}

J_2_kcal = 0.000239006
nm_2_Angstrom = 10.0


def get_value(symbol, key):
    return data.get(symbol, {}).get(key, None)

import numpy as np

import argparse
parser = argparse.ArgumentParser(description='generate MOMB pair style')
# Arguments supported by the code.
parser.add_argument('--atom_type', nargs='+')

args        = parser.parse_args()
atom_type       = args.atom_type

ntype = len(atom_type)

for i in range(ntype):
    for j in range(i,ntype):
        C6i = get_value(atom_type[i], "C6") *J_2_kcal * (nm_2_Angstrom ** 6)
        C6j = get_value(atom_type[j], "C6") *J_2_kcal * (nm_2_Angstrom ** 6)
        R0i = get_value(atom_type[i], "R0") 
        R0j = get_value(atom_type[j], "R0") 

        C6ij = np.sqrt(C6i*C6j)
        R0ij = R0i + R0j
        print("pair_coeff", i+1, j+1, "momb 0.0 1.0 1.0", f"{C6ij:.4f}", f"{R0ij:.4f}")


