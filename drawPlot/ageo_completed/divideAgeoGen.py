import numpy as np

# Load data from ageo.txt and gen.txt using NumPy
ageo_data = np.loadtxt('ageo_numi_right.txt')
gen_data = np.loadtxt('gen_numi_right.txt')

# Extract the first column from both ageo and gen
ageo_first_column = ageo_data[:, 0]
gen_first_column = gen_data[:, 0]

# Exclude the first column from both ageo and gen
ageo_data_exclude_first = ageo_data[:, 1:]
gen_data_exclude_first = gen_data[:, 1:]

# Perform division element-wise
result_data = ageo_data_exclude_first / gen_data_exclude_first

# Add the first column back to the result
result_data_with_first_column = np.column_stack((ageo_first_column, result_data))

# Save the result to a separate file
np.savetxt('ageo_numi_merged.txt', result_data_with_first_column, delimiter=' ')

