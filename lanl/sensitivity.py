import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib import rc

def generateSubmet(charge, mass):
    N_pot = 5 * 10 ** 21  # 0.4x0.4 m2 2 years
    acc_rate = 6 * 10 ** (-5)
    alpha = 1.0 / 137

    a = 2  # num of layers

    m_e = 0.00051
    m_pi = 0.135
    m_eta = 0.548
    m_rho = 0.775
    m_omega = 0.782
    m_phi = 1.019
    m_Jpsi = 3.1

    c_pi = 1.9
    c_eta = 0.21
    c_rho = 0.24
    c_omega = 0.24
    c_phi = 4.9e-03
    c_Jpsi = 5e-09

    branch_pi = 0.98
    branch_eta = 0.39
    branch_rho = 4.72e-5
    branch_omega = 7.28e-5
    branch_phi = 2.95e-4
    branch_Jpsi = 0.06

    N_gamma = 4.2 * 10 ** 5 * 0.3 * 5  # n_pe within pmt

    # N_gamma = 4.2*10**5*0.3*2 # n_pe within pmt, 20cm BGO
    # N_gamma = 4*10**5

    def I_2(x, y):
        return ((1 + 2 * x) * (1 - 4 * x) ** 0.5) / ((1 + 2 * y) * (1 - 4 * y) ** 0.5)

    def I_3(z, x):
        return 2 / (3 * 3.14) * ((1 - 4 * x / z) ** 0.5) * ((1 - z) ** 3) * (2 * x + z) / z ** 2

    def dl(x):
        return quad(I_3, 4 * x, 1, args=(x))[0]  # pass x to I
    v_dl = np.vectorize(dl)

    pi = np.zeros(np.shape(mass))
    eta = np.zeros(np.shape(mass))
    rho = np.zeros(np.shape(mass))
    omega = np.zeros(np.shape(mass))
    phi = np.zeros(np.shape(mass))
    jpsi = np.zeros(np.shape(mass))

    pi[mass < m_pi / 2] = (1 - np.exp(-N_gamma * charge[mass < m_pi / 2] ** 2)) ** a * 2 * c_pi * branch_pi * alpha * N_pot * acc_rate * v_dl(mass[mass < m_pi / 2] ** 2 / m_pi ** 2) * charge[mass < m_pi / 2] ** 2

    eta[mass < m_eta / 2] = (1 - np.exp(-N_gamma * charge[mass < m_eta / 2] ** 2)) ** a * 2 * c_eta * branch_eta * alpha * N_pot * acc_rate * v_dl(mass[mass < m_eta / 2] ** 2 / m_eta ** 2) * charge[mass < m_eta / 2] ** 2

    rho[mass < m_rho / 2] = (1 - np.exp(-N_gamma * charge[mass < m_rho / 2] ** 2)) ** a * 2 * c_rho * branch_rho * N_pot * acc_rate * I_2(mass[mass < m_rho / 2] ** 2 / m_rho ** 2, m_e ** 2 / m_rho ** 2) * charge[mass < m_rho/ 2] ** 2

    omega[mass < m_omega / 2] = (1 - np.exp(-N_gamma * charge[mass < m_omega / 2] ** 2)) ** a * 2 * c_omega * branch_omega * N_pot * acc_rate * I_2(mass[mass < m_omega / 2] ** 2 / m_omega ** 2, m_e ** 2 / m_omega ** 2) * charge[mass < m_omega / 2] ** 2

    phi[mass < m_phi / 2] = (1 - np.exp(-N_gamma * charge[mass < m_phi / 2] ** 2)) ** a * 2 * c_phi * branch_phi * N_pot * acc_rate * I_2(mass[mass < m_phi/ 2] ** 2 / m_phi ** 2, m_e ** 2 / m_phi ** 2) * charge[mass < m_phi / 2] ** 2

    jpsi[mass < m_Jpsi / 2] = (1 - np.exp(-N_gamma * charge[mass < m_Jpsi / 2] ** 2)) ** a * 2 * c_Jpsi * branch_Jpsi * N_pot * acc_rate * I_2(mass[mass < m_Jpsi / 2] ** 2 / m_Jpsi ** 2, m_e ** 2 / m_Jpsi ** 2) * charge[mass < m_Jpsi / 2] ** 2

    sensitivity = pi + eta + rho + omega + phi + jpsi

    return sensitivity

def generateLANL(charge, mass):
    N_pot = 5.9 * 10 ** 22 
    alpha = 1.0 / 137

    a = 2  # num of layers

    m_e = 0.00051
    m_pi = 0.135
    m_eta = 0.548

    c_pi = 1.9
    c_eta = 0.21
    c_rho = 0.24

    branch_pi = 0.98
    branch_eta = 0.39

    N_gamma = 4.2 * 10 ** 5 * 0.3 * 5  # n_pe within pmt

    ageo_data = np.loadtxt("ageo.txt")
    ageo10m = ageo_data[:, 1]
    ageo35m = ageo_data[:, 2]

    def I_3(z, x):
        return 2 / (3 * 3.14) * ((1 - 4 * x / z) ** 0.5) * ((1 - z) ** 3) * (2 * x + z) / z ** 2

    def dl(x):
        return quad(I_3, 4 * x, 1, args=(x))[0]  # pass x to I
    v_dl = np.vectorize(dl)

    pi = np.zeros(np.shape(mass))
    eta = np.zeros(np.shape(mass))

    pi[mass < m_pi / 2] = (1 - np.exp(-N_gamma * charge[mass < m_pi / 2] ** 2)) ** a * 2 * c_pi * branch_pi * alpha * N_pot * ageo10m[mass < m_pi/2] * v_dl(mass[mass < m_pi / 2] ** 2 / m_pi ** 2) * charge[mass < m_pi / 2] ** 2

    eta[mass < m_eta / 2] = (1 - np.exp(-N_gamma * charge[mass < m_eta / 2] ** 2)) ** a * 2 * c_eta * branch_eta * alpha * N_pot * ageo10m[mass < m_eta/2] * v_dl(mass[mass < m_eta / 2] ** 2 / m_eta ** 2) * charge[mass < m_eta / 2] ** 2
    
    sensitivity10m = pi + eta

    pi = np.zeros(np.shape(mass))
    eta = np.zeros(np.shape(mass))

    pi[mass < m_pi / 2] = (1 - np.exp(-N_gamma * charge[mass < m_pi / 2] ** 2)) ** a * 2 * c_pi * branch_pi * alpha * N_pot * ageo35m[mass < m_pi/2] * v_dl(mass[mass < m_pi / 2] ** 2 / m_pi ** 2) * charge[mass < m_pi / 2] ** 2

    eta[mass < m_eta / 2] = (1 - np.exp(-N_gamma * charge[mass < m_eta / 2] ** 2)) ** a * 2 * c_eta * branch_eta * alpha * N_pot * ageo35m[mass < m_eta/2] * v_dl(mass[mass < m_eta / 2] ** 2 / m_eta ** 2) * charge[mass < m_eta / 2] ** 2

    sensitivity35m = pi + eta

    return sensitivity10m, sensitivity35m

if __name__ == '__main__':
    mass   = np.logspace(-2, 1, 500)
    charge = np.logspace(-5, 1, 500)
    masses, charges = np.meshgrid(mass, charge)
    submet = generateSubmet(masses, charges)

    mass = np.loadtxt("mass.txt")
    masses, charges = np.meshgrid(mass, charges)
    lanl10m, lanl35m = generateLANL(masses, charges)

    

    fig, ax = plt.figure(figsize(16, 12), dpi=500)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]')
    plt.ylabel('$\epsilon=Q/e$')
    plt.xlim(0.01, 10)
    plt.ylim(0.00001, 1)

    plt.contour(mass, charge, submet, [12], colors = 'red', linestyles = 'dashed')



    plt.legend(loc=4)
    plt.legend(ncol=2)
# plt.savefig('sensitivity_lanl.pdf')
    plt.show()

