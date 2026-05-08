input_filename = "filtered_data.dat"
file_count = 1
current_chunk = []

with open(input_filename, "r") as infile:
    for line in infile:
        # Strip whitespace and check if the first column is 13.0000
        parts = line.split()
        if not parts:
            continue
            
        current_chunk.append(line)
        
        # Check if first column is 13.0000 (adjust decimals if needed)
        if float(parts[0]) == 13.0000:
            output_filename = f"split_{file_count}.txt"
            with open(output_filename, "w") as outfile:
                outfile.writelines(current_chunk)
            
            print(f"Created {output_filename} with {len(current_chunk)} lines.")
            file_count += 1
            current_chunk = []

# Optional: Handle any leftover data that didn't end with 13.0000
if current_chunk:
    with open(f"split_{file_count}_remnant.txt", "w") as outfile:
        outfile.writelines(current_chunk)
