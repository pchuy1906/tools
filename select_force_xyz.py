import sys

# sys.argv[0] is the script name. 
# We expect len(sys.argv) to be exactly 2 for "python script.py arg1"
if len(sys.argv) != 2:
    print("Error: Exactly one argument is required.", file=sys.stderr)
    print("Usage: python script.py <argument>", file=sys.stderr)
    sys.exit(1)

file_xyz = sys.argv[1]


import numpy as np

def filter_top_max_force_frames(file_path, top_n=20):
    frame_max_magnitudes = []
    
    with open(file_path, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            
            num_atoms = int(line.strip().split()[0])
            f.readline()  # Skip the header line
            
            frame_forces = []
            for _ in range(num_atoms):
                atom_line = f.readline().strip().split()
                # Parse out fx, fy, fz forces (columns 4, 5, 6)
                forces = [float(x) for x in atom_line[4:7]]
                frame_forces.append(forces)
            
            # Compute individual magnitudes within the frame
            magnitudes = np.linalg.norm(np.array(frame_forces), axis=1)
            
            # Select ONLY the largest force magnitude for this frame
            frame_max_magnitudes.append(np.max(magnitudes))
            
    # Convert to a stable NumPy array tracking 1 value per frame
    frame_max_magnitudes = np.array(frame_max_magnitudes)
    
    # Calculate the cutoff value for the top N frames
    # Partitioning is more performant than full sorting for large trajectories
    threshold = np.partition(frame_max_magnitudes, -top_n)[-top_n]
    
    # Create the logical (Boolean) array matching only your target frames
    logical_array = frame_max_magnitudes >= threshold
    
    # Extract only the max force values for the top 20 frames
    top_frame_values = frame_max_magnitudes[logical_array]
    top_frame_indices = np.where(logical_array)[0]
    
    return logical_array, top_frame_indices, top_frame_values

# Run analysis
logical_mask, indices, max_vals = filter_top_max_force_frames(file_xyz, top_n=20)

print(logical_mask)

def save_filtered_xyz(input_path, output_path, logical_mask):
    """
    Saves a new XYZ file containing only the frames where logical_mask is False.
    """
    frame_idx = 0

    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        while True:
            line = infile.readline()
            if not line:
                break  # End of file reached

            # Step 1: Read the atom count line
            num_atoms = int(line.strip().split()[0])

            # Step 2: Read the header line
            header_line = infile.readline()

            # Step 3: Read all atom lines for this specific frame
            atom_lines = [infile.readline() for _ in range(num_atoms)]

            # Step 4: Check the mask for the current frame
            # If False, write the complete frame chunk to the new file
            if not logical_mask[frame_idx]:
                outfile.write(line)          # Atom count line
                outfile.write(header_line)   # Header line
                outfile.writelines(atom_lines) # All atom data

            frame_idx += 1

out_file_xyz = f"filtered_{file_xyz}"
save_filtered_xyz(file_xyz, out_file_xyz, logical_mask)

