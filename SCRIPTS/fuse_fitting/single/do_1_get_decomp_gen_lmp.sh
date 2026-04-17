cat > aaa.py << EOF
import sys
from ase.io import read
from collections import defaultdict


lmp_data_file = sys.argv[1]
atoms = read(lmp_data_file, format='lammps-data', units='real')

v_mol_id = atoms.arrays["mol-id"]
v_type = atoms.arrays["type"]

# Grouping types by mol-id
mol_to_types = defaultdict(list)
for m_id, t in zip(v_mol_id, v_type):
    mol_to_types[m_id].append(t)
# To get just the list of lists:
final_list = list(mol_to_types.values())

max_mol = max(v_mol_id)
if len(final_list) != max_mol:
    raise ValueError(
        f"Molecule count mismatch! Expected {expected_count} molecules "
        f"(based on max mol-id), but found {actual_count} unique groups."
    )

string_vars = ["".join(map(str, mol)) for mol in final_list]
# print(string_vars)

count_first = string_vars.count(string_vars[0])
count_last = string_vars.count(string_vars[-1])

na, nb = count_first, count_last

if string_vars[0] == string_vars[-1]:
    na = count_last//2
    nb = count_last//2

with open("decomposition.dat", "w") as f:
    f.write(f"{na} {nb}\n")

EOF


cwd=`pwd`
for fold in $(ls -1vd LAMMPSSimulator/condensed_*_*_* ) ; do
    cd $fold
        pwd
        python $cwd/aaa.py  00000/LMP_condensed_data.lammps
        head decomposition.dat

    cd $cwd
done

