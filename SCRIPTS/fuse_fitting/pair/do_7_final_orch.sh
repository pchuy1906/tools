cat > script_orch.py << EOF
import os
import sys
import numpy as np
from ase.io import read
from orchestrator.utils.setup_input import setup_orch_modules
from orchestrator.utils.data_standard import (
    ENERGY_KEY,
    METADATA_PROPERTY_MAP,
    MOL_ID_KEY,
    MOL_TYPE_KEY,
    MOL_PROPERTY_MAP,
)


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
    'storage': {
        'storage_type': 'COLABFIT',
        'storage_args': {
            'credential_file':
            '/usr/workspace/fuse/orchestrator_storage_credential_fuse.json'
        },
    },
    'trainer': {
        "trainer_type": "Class1",
        "trainer_args": {
            "molecule_1": "tmp1",
            "molecule_2": "tmp2"
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


def save_dataset(selected_dimers, mol_name):
    """
    Save a dataset of dimer scan configurations to storage.
    """

    # Build base property map from metadata and molecular property definitions
    base_keys = {
        item["new_property_name"]: item["new_map"]
        for item in [METADATA_PROPERTY_MAP, MOL_PROPERTY_MAP]
    }

    # Add energy-specific mapping
    energy_keys = {"energy_field": ENERGY_KEY}

    # Set property maps in storage
    storage.set_property_map(keys=base_keys | energy_keys)

    # Generate a unique dataset name based on the dimer scan
    dataset_name = storage.generate_dataset_name(
        "dimer_scan",
        f"{len(selected_dimers)}",
        check_uniqueness=True,
    )

    # Assemble dataset metadata
    metadata = {
        "description": ("Dimer scan configurations generated for"
                        f"{mol_name}")
    }

    # Create the dataset and return its identifier
    dataset_id = storage.new_dataset(
        dataset_name,
        selected_dimers,
        metadata,
    )
    return dataset_id


def dimer_selection(dimers, de_dft, de_ff, energy_threshold=2000.0):
    """
    Select dimers based on a DFT binding energy threshold and compute
    the difference between DFT and force field binding energies.

    This function:
      1. Computes the energy difference: de = de_dft - de_ff.
      2. Filters dimers whose absolute DFT binding energy is below
         a specified threshold.
      3. Returns the selected dimers, their COM positions, and the
         energy differences for the selected set.

    Parameters
    ----------
    dimers : list
        List of dimer objects or configurations. The length must match
        the length of de_dft and de_ff.

    de_dft : array_like
        DFT binding energies for each dimer. Can be a list or NumPy array.

    de_ff : array_like
        Force field binding energies for each dimer. Same length as de_dft.

    energy_threshold : float, optional
        Absolute DFT energy cutoff used to select dimers. Only dimers
        with abs(de_dft[i]) < energy_threshold are kept.
        Default is 2000.0.

    Returns
    -------
    selected_dimers : list
        List of dimer objects that pass the energy threshold filter.

    e_diff_selected : np.ndarray
        1D NumPy array of energy differences (de_dft - de_ff) for the
        selected dimers.
    """
    print("\nDimer selection")
    print(f"  number of dimers before: {len(dimers)}")

    # Convert energies to NumPy arrays for safe vectorized operations
    de_dft = np.asarray(de_dft, dtype=float)
    de_ff = np.asarray(de_ff, dtype=float)

    # Basic shape consistency check to avoid subtle bugs
    if de_dft.shape != de_ff.shape:
        raise ValueError(
            f"Shape mismatch: de_dft {de_dft.shape} vs de_ff {de_ff.shape}")
    if len(dimers) != de_dft.shape[0]:
        raise ValueError(
            f"mismatch: dimers {len(dimers)} vs de_dft {de_dft.shape[0]}")

    # Energy difference between DFT and force field
    de = de_dft - de_ff

    selected_dimers = []
    e_diff_selected = []

    # Iterate over all dimers and select based on the DFT energy threshold
    for idx in range(len(de_dft)):
        if abs(de_dft[idx]) < energy_threshold:
            selected_dimers.append(dimers[idx])
            e_diff_selected.append(de[idx])

    print(f"  number of dimers after: {len(selected_dimers)}")

    # Convert lists to NumPy arrays for convenient downstream processing
    e_diff_selected = np.array(e_diff_selected, dtype=float)

    return selected_dimers, e_diff_selected


def gen_dimers(mol1, mol2):
    """
    Generate dimer configurations by rotating one molecule and varying
    the center-of-mass separation along a fixed direction.
    """
    # Rotation step in degrees
    angle_step_deg = 30
    rotation_angles = list(range(0, 360, angle_step_deg))
    n_rotations = len(rotation_angles)

    # rotations in degs - only around x axis
    z_rots = np.array(rotation_angles).reshape(-1, 1)
    mol1_rotations = np.hstack((np.zeros([z_rots.size, 2]), z_rots))
    mol2_rotations = np.array([0, 0, 180])

    # translation vectors - 1D vector, no iteration for now
    mol1_translation = np.zeros(3)
    mol2_translation = np.zeros(3)

    # distance between mol COMs in A - displacement along x dir
    radial_translation_dir = np.array([1, 0, 0])
    radial_displacements = np.linspace(3, 13, 51)

    # generate dimer structures
    dimers = augmentor.generate_dimer_configs(
        mol1,
        mol2,
        mol1_rotations,
        mol1_translation,
        mol2_rotations,
        mol2_translation,
        radial_translation_dir,
        radial_displacements,
    )

    # Center-of-mass distances, one for each dimer configuration
    rcoms = np.tile(radial_displacements, n_rotations)

    return dimers, rcoms

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
mol_name = f"{mol_a}_{mol_b}"

output_opt_ff = "/p/lustre5/pham20/workdir_1/orchestrator/examples/class_1_FF/obtained_opt_ff"

# Required input files for this pair
data_file_1 = f"{output_opt_ff}/opt_data_{mol_a}.lammps"
data_file_2 = f"{output_opt_ff}/opt_data_{mol_b}.lammps"
style_file_1 = f"{output_opt_ff}/style_{mol_a}.json"
style_file_2 = f"{output_opt_ff}/style_{mol_b}.json"

topo_mol_1 = potential._parse_to_topology(data_file_1, style_file_1)
topo_mol_2 = potential._parse_to_topology(data_file_2, style_file_2)

monomer_1 = potential.extract_monomer(topo_mol_1)
monomer_2 = potential.extract_monomer(topo_mol_2)
monomer_2.arrays[MOL_ID_KEY] = np.full(len(monomer_2), 2)
monomer_2.arrays[MOL_TYPE_KEY] += np.max(monomer_1.arrays[MOL_TYPE_KEY])

dimers, rcoms = gen_dimers(monomer_1, monomer_2)

path_de_lmp = f"LAMMPSSimulator/dimer_{mol_a}_{mol_b}/EINT.txt"
de_ff = np.loadtxt(path_de_lmp)

path_de_vasp = f"VaspOracle/dimer_{mol_a}_{mol_b}/EINT.txt"
de_dft = np.loadtxt(path_de_vasp)


selected_dimers, e_diff_dimer = dimer_selection(
    dimers,
    de_dft,
    de_ff,
    energy_threshold=2.0,
)

for dimer, energy in zip(selected_dimers, e_diff_dimer):
    dimer.info[ENERGY_KEY] = energy


print("\nGenerate and compute energies for condensed configurations")
all_condenseds = []
all_energies = []
mole_ratios = [0.1, 0.3, 0.5, 0.7, 0.9]
densities = [1000, 1200, 1500]
for mole_ratio in mole_ratios:
    for density in densities:
        condensed, decomposition = augmentor.generate_condensed_configs(
            monomer_1,
            monomer_2,
            mole_ratio=mole_ratio,
            cubic_length=20.0,
            density=density,
        )

        path_de_lmp = f"LAMMPSSimulator/condensed_{mol_a}_{mol_b}_{mole_ratio}_{density}/EINT.txt"
        de_ff = np.loadtxt(path_de_lmp)
        path_de_vasp = f"VaspOracle/condensed_{mol_a}_{mol_b}_{mole_ratio}_{density}/EINT.txt"
        de_dft = np.loadtxt(path_de_vasp)

        dump_file = f"LAMMPSSimulator/condensed_{mol_a}_{mol_b}_{mole_ratio}_{density}/00000/initial.dump"
        dump_frames = read(dump_file, format="lammps-dump-text", index=":")

        print(f"\nmolecule ratios {mole_ratio} ")
        print(f"and density {density}(kg/m3) ")
        print("Frames in dump:", len(dump_frames))

        all_configs = []
        for i, frame in enumerate(dump_frames):
            if len(frame) != len(condensed):
                raise ValueError(f"Frame {i}: atom count mismatch, "
                                 f"ref={len(condensed)}, frame={len(frame)}")
            atoms_i = condensed.copy()
            atoms_i.positions[:] = frame.positions
            atoms_i.cell = frame.cell
            atoms_i.pbc = frame.pbc
            if abs(de_ff[i]) < 1000.0:
                all_configs.append(atoms_i)
                all_energies.append(de_dft[i]-de_ff[i])

        print(len(all_configs))

        all_condenseds.extend(all_configs)



for condensed, energy in zip(all_condenseds, all_energies):
    condensed.info[ENERGY_KEY] = energy

print(len(all_condenseds), len(all_energies))

training = selected_dimers + all_condenseds
dataset_id = save_dataset(training, mol_name)
path_type = f"workdir_{mol_name}"

ij_pair_coeffs = trainer.train(
    path_type=path_type,
    potential=potential,
    storage=storage,
    dataset_list=[dataset_id],
    same_molecule=False,
)

ij_pair_file = f"opt_{mol_name}.json"
potential._write_ij_pair(
    ij_pair_file,
    ij_pair_coeffs,
    shift_mol_2=True,
)

EOF


fold_dimer=`ls -d LAMMPSSimulator/dimer_*`
echo $fold_dimer

python script_orch.py $fold_dimer

