import random

# Input and output file names
input_file = "POSCAR"
output_file = "POSCAR_perturbed"

# Read the POSCAR
with open(input_file, "r") as f:
    lines = f.readlines()

# Identify where coordinates start
# POSCAR format: coordinates start after 7th line (0-based index 7)
# But can vary if element line or number line differ; webll detect bDirectb or bCartesianb
for i, line in enumerate(lines):
    if line.strip().lower().startswith("direct") or line.strip().lower().startswith("cartesian"):
        coord_start = i + 1
        break

# Perturb the coordinates
perturbed_lines = lines[:coord_start]  # keep header unchanged

for line in lines[coord_start:]:
    if not line.strip():
        continue
    parts = line.split()
    if len(parts) < 3:
        perturbed_lines.append(line)
        continue
    x, y, z = map(float, parts[:3])
    # add random perturbation between -0.01 and 0.01
    x += random.uniform(-0.01, 0.01)
    y += random.uniform(-0.01, 0.01)
    z += random.uniform(-0.01, 0.01)
    perturbed_lines.append(f"  {x: .12f}  {y: .12f}  {z: .12f}\n")

# Write the new file
with open(output_file, "w") as f:
    f.writelines(perturbed_lines)

print(f"Perturbed POSCAR written to {output_file}")
