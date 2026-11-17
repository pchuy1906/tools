import sys


def count_chno_in_xyz(file_path):
    # Initialize counts for C, H, N, and O
    counts = {"C": 0, "H": 0, "N": 0, "O": 0}

    try:
        with open(file_path, "r") as f:
            lines = f.readlines()

        # Process atom data starting from line 2
        for line in lines[2:]:
            parts = line.strip().split()
            if parts:
                element = parts[0]
                if element in counts:
                    counts[element] += 1

        return counts

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        sys.exit(1)


if __name__ == "__main__":
    # Check if the user provided the filename argument
    if len(sys.argv) < 2:
        print("Usage: python main.py <filename.xyz>")
        sys.exit(1)

    # Get the filename from the command line argument
    filename = sys.argv[1]

    # Run the counter
    atom_counts = count_chno_in_xyz(filename)

    # Calculate total elements found
    total_atoms = sum(atom_counts.values())

    print(f"Atom Counts & Percentages for {filename}")
    if total_atoms > 0:
        for element, count in atom_counts.items():
            percentage = (count / total_atoms) * 100
            print(f"  {element}: {count} ({percentage:.0f}%)")
        print(f"  Total CHNO Atoms: {total_atoms}")
    else:
        print("  No C, H, N, or O atoms found in the file.")

