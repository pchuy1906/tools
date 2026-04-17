import sys
import numpy as np

input_filename = sys.argv[1]
output_filename = "filtered_data.dat"

data = []

# 1. Read and filter the data
try:
    with open(input_filename, "r") as f:
        for line in f:
            parts = line.split()
            
            # Process only lines with exactly 4 columns
            if len(parts) == 4 and abs(float(parts[1]))<=2.0:
                try:
                    # Convert each part to a float
                    row = [float(val) for val in parts]
                    data.append(row)
                except ValueError:
                    # Skip lines that have 4 words but aren't numbers (e.g. headers)
                    continue
except FileNotFoundError:
    raise FileNotFoundError(f"Could not find the file: {input_filename}")

# 2. Write the filtered data to the new file
with open(output_filename, "w") as f:
    for row in data:
        # Format: 12 spaces wide, 4 decimal places for alignment
        formatted_line = f"{row[0]:12.4f} {row[1]:12.4f} {row[2]:12.4f} {row[3]:12.4f}\n"
        f.write(formatted_line)

print(f"Processed {len(data)} lines. Cleaned data saved to '{output_filename}'.")

data_arr = np.array(data)

col2 = data_arr[:, 1]
col3 = data_arr[:, 2]
rmse = np.sqrt(np.mean((col2 - col3)**2))
print(f"RMSE between column 2 and 3: {rmse:.4f}")

col4 = data_arr[:, 3]
rmse = np.sqrt(np.mean((col2 - col4)**2))
print(f"RMSE between column 2 and 4: {rmse:.4f}")


