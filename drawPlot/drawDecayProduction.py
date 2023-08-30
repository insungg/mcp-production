import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter

from scipy.integrate import quad

# define matplotlib parameters
SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE + 5)  # fontsize of the figure title


def mesonDecayProduction(mass, NPOT):
	# EM constant
	alpha = 1.0 / 137

	# mass in GeV 
	m_e = 0.00051
    m_pi = 0.135
    m_eta = 0.548
	# m_rho = 0.775
	# m_omega = 0.782
	# m_phi = 1.019
    m_jpsi = 3.1
	m_upsilon = 9.46

	# meson / NPOT obtained from PYTHIA
    c_pi = 1.9
    c_eta = 0.21
	# c_rho = 0.24
	# c_omega = 0.24
	# c_phi = 4.9e-03
    c_jpsi = 3.81e-5
	c_upsilon = 2.5e-9

	# m -> e+e- branching ratio
    branch_pi = 0.98
    branch_eta = 0.39
	# branch_rho = 4.72e-5
	# branch_omega = 7.28e-5
	# branch_phi = 2.95e-4
    branch_jpsi = 0.05971
	branch_upsilon = 0.0238

	# read geometric acceptance obtained from PYTHIA
	ageo_pi = np.array([])
	ageo_eta = np.array([])
	ageo_jpsi = np.array([])
	ageo_upsilon = np.array([])
	with open('ageo.txt', 'r') as f:
		data = f.readlines()
		for line in data:
			values = line.strip().split()
			ageo_pi.append(float(values[1]))
			ageo_eta.append(float(values[2]))
			ageo_jpsi.append(float(values[3]))
			ageo_upsilon.append(float(values[4]))
			
			

	# define phase space integrals 
	def I2(x, y):
        return ((1 + 2 * x) * (1 - 4 * x) ** 0.5) / ((1 + 2 * y) * (1 - 4 * y) ** 0.5)

    def I3_integrand(z, x):
        return 2 / (3 * 3.14) * ((1 - 4 * x / z) ** 0.5) * ((1 - z) ** 3) * (2 * x + z) / (z ** 2)

    def I3(x):
        return quad(I3_integrand, 4 * x, 1, args=(x))[0]  # integrate
    v_I3 = np.vectorize(I3) # vectorize

	# define productions
	pi      = np.zeroes(np.shape(mass))
	eta     = np.zeroes(np.shape(mass))
	jpsi    = np.zeroes(np.shape(mass))
	upsilon = np.zeroes(np.shape(mass))

	# masking to achieve kinematical validity
	# Dalitz decay
	pi[ mass < m_pi/2 ] = NPOT * ageo_pi[mass < m_pi/2 ] * 2 * c_pi * branch_pi * alpha * v_I3( mass[ mass < m_pi/2 ] ** 2 / mass[ mass < m_pi/2 ] ** 2)
	eta[ mass < m_eta/2 ] = NPOT * ageo_eta[mass < m_eta/2 ] * 2 * c_eta * branch_eta * alpha * v_I3( mass[ mass < m_eta/2 ] ** 2 / mass[ mass < m_eta/2 ] ** 2)
	# Direct decay
	jpsi[ mass < m_jpsi/2 ] = NPOT * ageo_jpsi[mass < m_jpsi/2 ] * 2  * c_jpsi * branch_jpsi * I2( mass[ mass < m_jpsi < 2] ** 2 / m_jpsi ** 2, m_e ** 2 / m_jpsi ** 2  )
	upsilon[ mass < m_upsilon/2 ] = NPOT * ageo_upsilon[mass < m_upsilon/2] * 2 * c_upsilon * branch_upsilon * I2( mass[ mass < m_upsilon < 2] ** 2 / m_upsilon ** 2, m_e ** 2 / m_upsilon ** 2  )
 
	production = pi + eta + jpsi + upsilon

	return production





if __name__ == '__main__':
	mass = np.array([])

	with open('mass.txt', 'r') as f:
		data = f.readline()
		for line in data:
			values = line.strip().split()
			mass.append(flaot(values[0]))

	numi = 6e20
	dune = 1e21
	chi_numi = mesonDecayProduction(mass, numi)
	chi_dune = mesonDecayProduction(mass, dune)

	fig = plt.figure( figsize=(13, 26))
	plt.xscale('log')
	plt.yscale('log')


	
