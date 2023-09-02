import numpy as np

a = np.logspace(-2, 1, 30)

np.savetxt("mass.txt", a)
