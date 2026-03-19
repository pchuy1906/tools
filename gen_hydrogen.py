import numpy as np
import sys

def read_xyz_from_selected_atoms(filename):
    xyz = []  # list of (element, x, y, z)

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if not line.startswith("Selected Atom"):
                continue

            parts = line.split()
            # Example tokens:
            # ["Selected", "Atom", "No.1:", "122", "O", "8", "+3.254020000", "+15.973611000", "+8.806910000"]

            element = parts[4]
            x = float(parts[6])
            y = float(parts[7])
            z = float(parts[8])

            xyz.append((x, y, z))

    return np.array(xyz)


def point_at_distance_np(xyz1, xyz2, d=1.0):
    """
    xyz1, xyz2: NumPy arrays of shape (3,) for points A and B.
    Returns H such that:
        - H lies on the line from A to B
        - distance(A, H) = d
    """
    xyz1 = np.asarray(xyz1, dtype=float)
    xyz2 = np.asarray(xyz2, dtype=float)

    v = xyz2 - xyz1              # vector AB
    L = np.linalg.norm(v)        # length of AB
    if L == 0:
        raise ValueError("xyz1 and xyz2 are identical, direction undefined")

    u = v / L                    # unit vector from A to B
    H = xyz1 + d * u             # point at distance d from A
    return H


if len(sys.argv) < 3:
    print("Usage: python script.py  input_file  distOH")
    sys.exit()

input_file = sys.argv[1]
distOH = sys.argv[2]
distOH = float(distOH)

xyz = read_xyz_from_selected_atoms(input_file)
nxyz = len(xyz)
nhydrogen = nxyz//2
for i in range(nhydrogen):
    xyz1 = xyz[2*i,:]
    xyz2 = xyz[2*i+1,:]
    H = point_at_distance_np(xyz1, xyz2, d=distOH)
    print(f"{'H':>5} " + " ".join(f"{x:12.6f}" for x in H))


