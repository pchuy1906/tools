import numpy as np
import sys

# Atomic masses (extend as needed)
atomic_masses = {
    'H': 1.008,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    # Add more elements here
}

def read_xyz_frames(filename):
    frames = []
    energies = []
    with open(filename) as f:
        lines = f.readlines()
    idx = 0
    while idx < len(lines):
        num_atoms = int(lines[idx].strip())
        atoms = []
        parts = lines[idx + 1].split()
        energies.append( float(parts[-1]) )
        for i in range(num_atoms):
            parts = lines[idx + 2 + i].split()
            symbol = parts[0]
            coords = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            if int(parts[-1])==2:
                atoms.append((symbol, coords))
        frames.append(atoms)
        idx += 2 + num_atoms
    return frames, energies

def center_of_mass(atoms):
    total_mass = 0.0
    com = np.zeros(3)
    for symbol, pos in atoms:
        mass = atomic_masses.get(symbol, 0.0)
        com += mass * pos
        total_mass += mass
    if total_mass == 0.0:
        return np.array([np.nan, np.nan, np.nan])
    return com / total_mass

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py file.xyz")
        sys.exit(1)
    xyz_file = sys.argv[1]
    frames, energies = read_xyz_frames(xyz_file)
    energies = np.array(energies)
    energies = energies - energies[-1]

    output_filename = "output.txt"
    with open(output_filename, "w") as f:
        for i, frame in enumerate(frames):
            com = center_of_mass(frame)
            f.write(f"{com[0]:12.4f} {energies[i]:12.4f}\n")
