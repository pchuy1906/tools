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


cat > bbb.py << EOF
import os
import sys
import numpy as np
from orchestrator.utils.setup_input import setup_orch_modules


input_args_orch = {
    'augmentor': {
        'augmentor_type': 'BASE',
        'augmentor_args': {},
    },
    'potential': {
        "potential_type": "Class1",
        "potential_args": {
            "kim_api": "kim-api-collections-management",
            "model_driver": "no-driver",
            "n_type_max": 20,
            "rcut": 10.0,
            "species": ["C", "H", "N", "O"]
        }
    },
}

(
    augmentor,
    descriptor,
    oracle,
    potential,
    score,
    simulator,
    storage,
    target_property,
    trainer,
    workflow,
) = setup_orch_modules(input_args_orch)

def move_files_lammps(lammps_files, pre_ind):
    """
    Rename LAMMPS files by appending a suffix using Unix mv via os.system.

    Parameters
    ----------
    lammps_files : iterable of str
        List or other iterable of file paths to rename.
    pre_ind : str or int
        Suffix to append to each filename.
    """
    for file in lammps_files:
        cmd = f"mv {file} {pre_ind}_{file}"
        status = os.system(cmd)
        if status != 0:
            print(f"Warning: failed to move {file} to {pre_ind}_{file}")

data_file = "data.lammps"
style_file = "style.lammps"
pair_file = "pair.lammps"
lammps_files = [data_file, style_file, pair_file]

path = sys.argv[1]

filename = path.split('/')[-1]
parts = filename.split('_')
mol_a = parts[1]
mol_b = parts[2]

output_opt_ff = "/p/lustre5/pham20/workdir_1/orchestrator/examples/class_1_FF/obtained_opt_ff"

# Required input files for this pair
data_file_1 = f"{output_opt_ff}/opt_data_{mol_a}.lammps"
data_file_2 = f"{output_opt_ff}/opt_data_{mol_b}.lammps"
style_file_1 = f"{output_opt_ff}/style_{mol_a}.json"
style_file_2 = f"{output_opt_ff}/style_{mol_b}.json"

topo_mol_1 = potential._parse_to_topology(data_file_1, style_file_1)
topo_mol_2 = potential._parse_to_topology(data_file_2, style_file_2)

# Write and relocate LAMMPS files for monomer 1
potential._write_potential_to_file(
    topo_mol_1,
    '.',
    write_pair_coeffs=False,
    write_ij_pair_coeffs=False,
)
move_files_lammps(lammps_files, pre_ind="LMP_monomer_1")

# Write and relocate LAMMPS files for monomer 2
potential._write_potential_to_file(
    topo_mol_2,
    '.',
    write_pair_coeffs=False,
    write_ij_pair_coeffs=False,
)
move_files_lammps(lammps_files, pre_ind="LMP_monomer_2")

multiplicity = np.loadtxt("decomposition.dat", dtype=int)
lj_inter_param_type="zeros"
# Merge monomer force fields using the specified LJ interaction parameters
merged_ff = potential.merge_ff(
    [topo_mol_1, topo_mol_2],
    lj_inter_param_type=lj_inter_param_type,
)

# Write and relocate LAMMPS files for the merged dimer force field
potential.write_topology_datafile(
    merged_ff,
    '.',
    multiplicity=multiplicity,
    write_pair_coeffs=False,
    write_ij_pair_coeffs=False,
)
if multiplicity is None:
    move_files_lammps(lammps_files, pre_ind="LMP_dimer")
else:
    move_files_lammps(lammps_files, pre_ind="LMP_condensed")


EOF



cwd=`pwd`
for fold in $(ls -1vd LAMMPSSimulator/condensed_*_*_*_* ) ; do
    cd $fold
        pwd
        python $cwd/aaa.py  00000/LMP_condensed_data.lammps
        python $cwd/bbb.py  $fold

    cd $cwd
done

