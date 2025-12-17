import numpy as np
import argparse
import sys

HA_TO_KCAL_MOL = 627.503
BOHR_TO_ANG = 0.529177

def parse_arguments():
    parser = argparse.ArgumentParser(description='Setup DFTB input for file_xyz')
    parser.add_argument("--file_xyz_1", default='file1.xyz', help='Input xyz file 1')
    parser.add_argument("--file_xyz_2", default='file2.xyz', help='Input xyz file 2')
    parser.add_argument("--cell_option_input", default='cell_3', help='Cell option input: "NON_ORTHO", "cell_3", or "cell_9"')
    parser.add_argument("--cell_option_output", default='cell_3', help='Cell option output: "NON_ORTHO", "cell_3", or "cell_9"')
    parser.add_argument("--export_quantities", default='xyzfes', help='Export quantities: "xyzfe" or "xyzfes"')
    return parser.parse_args()

def process_frames(f1, f2, f3, f4, f5, f6, cell_option_input, cell_option_output, export_quantities):
    nframe = 0
    while True:
        line1 = f1.readline()
        if not line1:
            break  # End of file

        natom = int(line1.strip())
        f2.readline()  # Skip corresponding line in file2
        f3.write(f"{natom}\n")

        cell1_data = f1.readline().split()
        cell2_data = f2.readline().split()

        if cell_option_input == "NON_ORTHO":
            if export_quantities == "xyzfes":
                cell1_9 = [float(x) for x in cell1_data[1:10]]
                stress1_6 = [float(x) for x in cell1_data[10:16]]
                energy1 = float(cell1_data[16])
                stress2_6 = [float(x) for x in cell2_data[10:16]]
                energy2 = float(cell2_data[16])
            elif export_quantities == "xyzfe":
                cell1_9 = [float(x) for x in cell1_data[1:10]]
                energy1 = float(cell1_data[10])
                energy2 = float(cell2_data[10])
                stress1_6 = [0.0] * 6
                stress2_6 = [0.0] * 6
            else:
                sys.exit("Export quantity not implemented for NON_ORTHO.")
        else:
            sys.exit(f"Cell option input '{cell_option_input}' not implemented.")

        dstress = np.subtract(stress1_6, stress2_6)
        denergy = energy1 - energy2

        if cell_option_output == "NON_ORTHO":
            f3.write("NON_ORTHO")
            for value in cell1_9:
                f3.write(f"{value:12.6f}")
            if export_quantities == "xyzfes":
                for value in dstress:
                    f3.write(f"{value:12.6f}")
            f3.write(f"{denergy:15.6f}\n")
        elif cell_option_output == "xyzfe":
            f3.write("NON_ORTHO")
            for value in cell1_9:
                f3.write(f"{value:12.6f}")
            f3.write(f"{denergy:15.6f}\n")
        else:
            sys.exit(f"Cell option output '{cell_option_output}' not implemented.")

        for k in range(natom):
            atom1_data = f1.readline().split()
            atom2_data = f2.readline().split()
            atom_symbol = atom1_data[0]
            xyz1 = [float(x) for x in atom1_data[1:4]]
            fxyz1 = [float(x) for x in atom1_data[4:7]]
            fxyz2 = [float(x) for x in atom2_data[4:7]]
            itype = int(atom1_data[7])
            imol = int(atom1_data[8])
            dfxyz = np.subtract(fxyz1, fxyz2)
            f3.write(
                f"{atom_symbol} {xyz1[0]:12.6f} {xyz1[1]:12.6f} {xyz1[2]:12.6f} "
                f"{dfxyz[0]:12.6f} {dfxyz[1]:12.6f} {dfxyz[2]:12.6f} {itype:5d} {imol:5d}\n"
            )
            for ixyz in range(3):
                f4.write(f"{fxyz1[ixyz] * HA_TO_KCAL_MOL / BOHR_TO_ANG:12.6f} FORCE_{nframe}_{atom_symbol}\n")
                f5.write(f"{fxyz2[ixyz] * HA_TO_KCAL_MOL / BOHR_TO_ANG:12.6f} FORCE_{nframe}_{atom_symbol}\n")

        f6.write(f"{natom}\n")
        f4.write(f"{energy1:12.6f} ENERGY_{nframe}\n")
        f5.write(f"{energy2:12.6f} ENERGY_{nframe}\n")

        nframe += 1

def main():
    args = parse_arguments()
    try:
        with open(args.file_xyz_1, "rt") as f1, \
             open(args.file_xyz_2, "rt") as f2, \
             open("ndelta.xyz", "w") as f3, \
             open("data_file1.dat", "w") as f4, \
             open("data_file2.dat", "w") as f5, \
             open("data_natom.dat", "w") as f6:
            process_frames(
                f1, f2, f3, f4, f5, f6,
                args.cell_option_input,
                args.cell_option_output,
                args.export_quantities
            )
    except FileNotFoundError as e:
        print(f"File error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
