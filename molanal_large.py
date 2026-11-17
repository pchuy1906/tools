import re
import numpy as np
from collections import defaultdict
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

def parse_bonds_dat(filename):
    print(f"Reading cutoffs from {filename}...")
    with open(filename, 'r') as f:
        content = f.read()
    
    # Matches any pair of letters followed by a float (e.g., H H 0.900 or C H 1.400)
    pattern = r"([A-Za-z]+)\s+([A-Za-z]+)\s+(\d+\.\d+)"
    matches = re.findall(pattern, content)
    
    bond_cutoffs = {}
    for el1, el2, cutoff_str in matches:
        cutoff_val = float(cutoff_str)
        bond_cutoffs[(el1, el2)] = cutoff_val
        
    print(f"Successfully parsed {len(bond_cutoffs)} bond cutoffs from {filename}")
    for pairs, val in bond_cutoffs.items():
        print(f"  {pairs[0]}-{pairs[1]}: {val} ")
        
    return bond_cutoffs

def split_and_unwrap_xyz_strict_cutoff(filename, bond_cutoffs):
    print(f"\nReading XYZ file: {filename}...")
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    num_atoms = int(lines[0].strip())
    comment_line = lines[1].strip()
    
    # Parse cell lengths (a, b, c) from the comment line
    cell_lengths = [float(val) for val in re.findall(r"[-+]?\d*\.\d+|\d+", comment_line)]
    if len(cell_lengths) < 3:
        raise ValueError(f"Could not parse 3 cell lengths from comment line: '{comment_line}'")
    
    box_size = np.array(cell_lengths[:3])
    print(f"Detected periodic box dimensions: {box_size}")
    
    # Parse elements and coordinates
    elements = []
    coords = np.zeros((num_atoms, 3))
    
    for idx, line in enumerate(lines[2:2+num_atoms]):
        parts = line.split()
        if len(parts) >= 4:
            elements.append(parts[0])  # Extract element symbol cleanly
            coords[idx] = [float(parts[1]), float(parts[2]), float(parts[3])]
            
    print("Wrapping raw coordinates into primary periodic box...")
    coords = coords % box_size
    
    print("Building periodic spatial tree...")
    # The maximum distance to query is strictly defined by the largest cutoff in bonds.dat
    if not bond_cutoffs:
        raise ValueError("No bond cutoffs found in bonds.dat. Cannot process molecules.")
        
    max_cutoff = max(bond_cutoffs.values())
    tree = cKDTree(coords, boxsize=box_size)
    candidate_pairs = tree.query_pairs(r=max_cutoff, output_type='ndarray')
    
    print(f"Evaluating {len(candidate_pairs)} candidate pairs with strict cutoffs...")
    if len(candidate_pairs) == 0:
        print("No bonds found.")
        return

    # Calculate exact periodic distances for candidate pairs
    u = candidate_pairs[:, 0]
    v = candidate_pairs[:, 1]
    diff = coords[u] - coords[v]
    diff = diff - box_size * np.round(diff / box_size)  # Minimum image convention
    dists = np.linalg.norm(diff, axis=1)
    
    # Vectorized element lookup
    elements_arr = np.array(elements)
    el_u = elements_arr[u]
    el_v = elements_arr[v]
    
    # Filter pairs using ONLY the defined pair cutoffs (unlisted pairs remain False)
    valid_mask = np.zeros(len(candidate_pairs), dtype=bool)
    
    for (e1, e2), cutoff_val in bond_cutoffs.items():
        # Match both orderings (e.g., C-H and H-C)
        mask = ((el_u == e1) & (el_v == e2)) | ((el_u == e2) & (el_v == e1))
        valid_mask[mask] = dists[mask] <= cutoff_val
        
    # Final filtered bonds list
    pairs = candidate_pairs[valid_mask]
    print(f"Found {len(pairs)} valid bonds strictly matching bonds.dat thresholds.")
    
    print("Identifying molecular networks...")
    row = pairs[:, 0]
    col = pairs[:, 1]
    data = np.ones(len(pairs))
    
    graph = csr_matrix((data, (row, col)), shape=(num_atoms, num_atoms))
    num_molecules, labels = connected_components(csgraph=graph, directed=False)
    print(f"Identified {num_molecules} distinct molecules.")
    
    # Group atom indices by molecule
    molecule_groups = [[] for _ in range(num_molecules)]
    for atom_idx, mol_label in enumerate(labels):
        molecule_groups[mol_label].append(atom_idx)
    
    # Build an adjacency list for unwrapping
    adj_list = defaultdict(list)
    for u_b, v_b in pairs:
        adj_list[u_b].append(v_b)
        adj_list[v_b].append(u_b)
        
    print("Unwrapping molecules across periodic boundaries for VMD...")
    unwrapped_coords = np.copy(coords)
    visited_unwrap = np.zeros(num_atoms, dtype=bool)
    
    for indices in molecule_groups:
        if not indices:
            continue
        root = indices[0]
        queue = [root]
        visited_unwrap[root] = True
        
        while queue:
            curr = queue.pop(0)
            curr_pos = unwrapped_coords[curr]
            
            for neighbor in adj_list[curr]:
                if not visited_unwrap[neighbor]:
                    diff_unwrapped = coords[neighbor] - curr_pos
                    diff_unwrapped = diff_unwrapped - box_size * np.round(diff_unwrapped / box_size)
                    unwrapped_coords[neighbor] = curr_pos + diff_unwrapped
                    visited_unwrap[neighbor] = True
                    queue.append(neighbor)
                    
    print("Writing unwrapped molecule files...")
    for idx, atom_indices in enumerate(molecule_groups):
        out_name = f'molecule_{idx+1}.xyz'
        with open(out_name, 'w') as f:
            f.write(f"{len(atom_indices)}\n")
            f.write(f"Molecule {idx+1} (Strict cutoffs from bonds.dat) | Box: {comment_line}\n")
            for i in atom_indices:
                pos = unwrapped_coords[i]
                f.write(f"{elements[i]} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
                
    print("Done! All strictly-defined molecules extracted and unwrapped.")

# ==========================================
# Execution Engine
# ==========================================
try:
    # 1. Parse the bonds configuration data file
    custom_cutoffs = parse_bonds_dat('bonds.dat')
    
    # 2. Run the strict molecule split process
    split_and_unwrap_xyz_strict_cutoff('last.xyz', bond_cutoffs=custom_cutoffs)

except Exception as e:
    print(f"\nAn error occurred: {e}")

