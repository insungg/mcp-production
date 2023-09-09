import numpy as np

# Load data from a.txt and b.txt into NumPy arrays
data_a = np.loadtxt('ageo_dy_dune.txt')
data_b = np.loadtxt('ageo_dy_dune_heavy.txt')

# Append data_b to data_a
combined_data = np.vstack((data_a, data_b))

# Sort the combined data by the first column (mass)
sorted_data = combined_data[combined_data[:, 0].argsort()]

# Save the sorted data to a new file
np.savetxt('ageo_dy_dune_full.txt', sorted_data)

# Print the sorted data
print(sorted_data)

