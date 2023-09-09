# Open the input file for reading
with open('mass.txt', 'r') as input_file:
    # Open a new file for writing the output
    with open('output.txt', 'w') as output_file:
        # Read each line from the input file along with its index
        for index, line in enumerate(input_file, start=1):
            # Strip any leading/trailing whitespace
            line = line.strip()
            
            # Check if the line is not empty
            if line:
                try:
                    # Attempt to convert the line to a floating-point number (mass)
                    mass = float(line)
                    # Write the output line with the formatted mass value and index
                    line_launch = f'launch -n mass_{index}\n'
                    output_file.write(line_launch)
                    line_setmass = f'set mchi {mass}\n'
                    output_file.write(line_setmass)
                except ValueError:
                    # Handle the case where the line is not a valid number
                    print(f"Skipping invalid line: {line}")

print("Conversion complete. Output saved to 'output.txt'.")

