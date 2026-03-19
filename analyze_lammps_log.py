import sys
import numpy as np

file_lmp_log="log.lammps"


def find_units(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line_number, line in enumerate(f, start=1):
                # case-sensitive: "units" exactly as typed
                if "units" in line:
                    print(f"\nLAMMPS {line.rstrip()}")
                    return line.split()[1]
        print("No line containing 'units' was found.")
    except FileNotFoundError:
        print(f"File not found: {path}")
    except Exception as e:
        print(f"An error occurred: {e}")


def find_printed_quantities(path):
    last_line = None
    last_line_number = None

    try:
        with open(path, "r", encoding="utf-8") as f:
            for line_number, line in enumerate(f, start=1):
                if "Step" in line:  # case-sensitive
                    last_line = line
                    last_line_number = line_number

        if last_line is not None:
            print("\nprinted quantities:")
            print(f"{last_line.rstrip()}")
            return last_line_number + 1, last_line.split()
        else:
            print("No line containing 'Step' was found.")
    except FileNotFoundError:
        print(f"File not found: {path}")
    except Exception as e:
        print(f"An error occurred: {e}")


def find_line_end(path):
    last_line_number = None
    last_line = None

    try:
        with open(path, "r", encoding="utf-8") as f:
            for line_number, line in enumerate(f, start=1):
                last_line_number = line_number
                last_line = line
                if "Loop" in line:  # case-sensitive
                    loop_line_number = line_number
                    loop_line = line

        if last_line_number is None:
            # Empty file
            print("File is empty.")
            return None

        if "loop_line_number" in locals():
            print("\nline end::")
            print(f"{loop_line.rstrip()}")
            return loop_line_number - 1
        else:
            # No 'Loop' found, return last line number
            print("No line containing 'Loop' was found.")
            print("Returning the last line number of the file.")
            print(f"Last line: {last_line.rstrip()}")
            return last_line_number-2

    except FileNotFoundError:
        print(f"File not found: {path}")
    except Exception as e:
        print(f"An error occurred: {e}")


def load_block_to_array(file_path, line_start, line_end):
    """
    Read lines [line_start, line_end] inclusive (1-based) from file_path
    and return them as a 2D NumPy array of floats.
    """
    data_rows = []
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            for i, line in enumerate(f, start=1):
                if i < line_start:
                    continue
                if i > line_end:
                    break
                stripped = line.strip()
                if not stripped:
                    # skip empty lines
                    continue
                # split on any whitespace and convert to float
                parts = stripped.split()
                row = [float(x) for x in parts]
                data_rows.append(row)

        if not data_rows:
            print("No data lines found in the specified range.")
            return None

        arr = np.array(data_rows, dtype=float)
        print("Shape of array:", arr.shape)
        return arr

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except ValueError as e:
        print("Conversion to float failed, check that the lines contain numeric data only.")
        print("Error:", e)
    except Exception as e:
        print("Unexpected error:", e)

lmp_units = find_units(file_lmp_log)

line_start, printed_quantities = find_printed_quantities(file_lmp_log)

line_end = find_line_end(file_lmp_log)


print (line_start, line_end)

data = load_block_to_array(file_lmp_log, line_start, line_end)

x_quantity = sys.argv[1]

print(f"x_axis is {x_quantity}")

if x_quantity in printed_quantities:
    idx = printed_quantities.index(x_quantity)
else:
    raise ValueError(f"{x_quantity} not found in printed_quantity")

f_convert_press = 1.0
if lmp_units=="real":
    f_convert_press = 0.000101325
keys_press = ["Px", "Py", "Pz", "Press"]

for i in range(len(printed_quantities)):
    quantity = printed_quantities[i]
    file_name = f"AAA_{quantity}.dat"

    f_convert = 1.0

    has_press = any(t in quantity for t in keys_press)
    if has_press:
        f_convert = f_convert_press

    y = data[:,i] * f_convert
    print(file_name)

    with open(file_name, "w") as f:
        for a, b in zip(x, y):
            f.write(f"{a} {b}\n")



