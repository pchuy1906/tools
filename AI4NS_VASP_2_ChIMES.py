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

def read_poscar(filename):
    """
    Reads atomic species and their counts from a VASP POSCAR file.

    Parameters:
        filename (str): Path to the POSCAR file.

    Returns:
        tuple: (species, counts)
            species (list of str): Atomic species.
            counts (list of int): Number of atoms for each species.
    """
    with open(filename, 'rt') as file:
        # Skip the first five lines (header information)
        for _ in range(5):
            file.readline()
        # Read atomic species
        species = file.readline().split()
        ntype = len(species)
        # Read atomic counts
        counts = [int(x) for x in file.readline().split()[:ntype]]
    return species, counts

def create_atom_list(species, counts):
    """
    Generates a list of atomic species labels, repeated according to their counts.

    Parameters:
        species (list of str): List of atomic species labels.
        counts (list of int): Corresponding number of atoms for each species.

    Returns:
        list of str: Expanded list of atomic species labels.
    """
    atom_list = []
    for label, count in zip(species, counts):
        atom_list.extend([label] * count)
    return atom_list

def read_outcar(filename, natom):
    """
    Extracts cell vectors, atomic positions, forces, energy, and stress tensor from a VASP OUTCAR file.

    Parameters:
        filename (str): Path to the OUTCAR file.
        natom (int): Number of atoms in the system.

    Returns:
        tuple:
            cell_vectors (np.ndarray): Cell lattice vectors (shape: (3, 3)).
            positions (np.ndarray): Atomic positions (shape: (natom, 3)).
            forces (np.ndarray): Atomic forces (shape: (natom, 3)), converted to Hartree/Bohr.
            energy (float): Free energy in kcal/mol.
            stresses (list of float): Stress tensor components in GPa.
    """
    positions = np.zeros((natom, 3))
    forces = np.zeros((natom, 3))
    cell_vectors = np.zeros((3, 3))
    energy = None
    stresses = None

    with open(filename, 'rt') as file:
        for line in file:
            # Extract stress tensor
            if "in kB" in line:
                tokens = line.split()
                stresses = [float(tokens[i]) * convert_kbar_2_GPa for i in range(2, 8)]
            # Extract cell lattice vectors
            elif "direct lattice vectors" in line:
                for i in range(3):
                    vector_line = next(file).split()
                    cell_vectors[i] = [float(vector_line[j]) for j in range(3)]
            # Extract atomic positions and forces
            elif "POSITION" in line:
                next(file)  # Skip header line
                for i in range(natom):
                    atom_line = next(file).split()
                    positions[i] = [float(atom_line[j]) for j in range(3)]
                    forces[i] = [float(atom_line[j]) for j in range(3, 6)]
                # Convert forces to Hartree/Bohr
                forces *= convert_eV_2_Hartree / convert_Angstrom_2_Bohr
            # Extract free energy
            elif "free  energy" in line:
                energy = float(line.split()[4]) * convert_eV_2_kcal_per_mol

    return cell_vectors, positions, forces, energy, stresses

def write_chimes_xyzf(
    atom_list,
    cell_type,
    cell_vectors,
    export_stress,
    energy,
    positions,
    forces,
    stresses=None,
    filename="VASP_2_ChIMES.xyzf",
    lammps_ids=None
):
    """
    Writes atomic configuration, cell parameters, energy, forces, and optional stress to a ChIMES-format .xyzf file.

    Parameters:
        atom_list (list of str): List of atomic species labels.
        cell_type (str): Cell type ('cell_3' or 'NON_ORTHO').
        cell_vectors (np.ndarray): Cell lattice vectors (shape: (3, 3)).
        export_stress (bool): Whether to include stress tensor in the output.
        energy (float): Total energy value.
        positions (np.ndarray): Atomic positions (shape: (n_atoms, 3)).
        forces (np.ndarray): Atomic forces (shape: (n_atoms, 3)).
        stresses (list of float, optional): Stress tensor components (length 6).
        filename (str): Output filename (default: "VASP_2_ChIMES.xyzf").
    """
    natom = len(atom_list)

    if lammps_ids is not None:
        sorted_data = lammps_ids[lammps_ids[:, 1].argsort()]
        order = sorted_data[:,0]
        atype = sorted_data[:,2]
        amol  = sorted_data[:,3]


    with open(filename, "w") as file:
        # Write number of atoms
        file.write(f"{natom}\n")

        # Write cell parameters
        if cell_type == "cell_3":
            for i in range(3):
                for j in range(3):
                    file.write(f"{cell_vectors[i, j]:15.9f}")
        elif cell_type == "NON_ORTHO":
            file.write("NON_ORTHO ")
            for i in range(3):
                for j in range(3):
                    file.write(f"{cell_vectors[i, j]:15.9f}")

        # Write stresses if requested
        if export_stress:
            if stresses is None or len(stresses) != 6:
                raise ValueError("Stresses must be provided as a list of 6 values when export_stress is True.")
            for value in stresses:
                file.write(f"{value:15.9f}")

        # Write energy
        file.write(f"{energy:20.9f}\n")

        # Write atomic data
        for iatom in range(natom):

            if lammps_ids is not None:
                i = order[iatom]
            else:
                i = iatom

            file.write(f"{atom_list[i]}")
            for j in range(3):
                file.write(f"{positions[i, j]:15.9f}")
            for j in range(3):
                file.write(f"{forces[i, j]:15.9f}")

            if lammps_ids is not None:
                file.write(f" {atype[iatom]} {amol[iatom]}")
            file.write("\n")

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

