from ase.io import read
import numpy as np

# Read all configurations
images = read('XDATCAR', index=':')

n_total = len(images)
n_select = 10

# Get 10 equally spaced indices
indices = np.linspace(0, n_total - 1, n_select, dtype=int)

# Print the indices of the selected configurations
print(f"Number of configurations: {len(images)}")
print("Selected configuration indices (order numbers):", indices.tolist())

with open('selected.xyz', 'w') as f:
    for i in indices:
        atoms = images[i]
        # Number of atoms
        f.write(f"{len(atoms)}\n")
        # Lattice parameters as comment line
        cell = atoms.cell.array.reshape(-1)
        f.write(' '.join(f"{x:.10f}" for x in cell) + '\n')
        # Atomic symbols and positions
        for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.positions):
            f.write(f"{symbol} {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}\n")
