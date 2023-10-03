import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib import rc
from matplotlib.ticker import MaxNLocator

import matplotlib
# matplotlib.use('TkAgg')

def generateSubmet(mass, charge):
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

    N_gamma = 2.5e5

    #N_gamma = 4.2 * 10 ** 5 * 0.3 * 5  # n_pe within pmt
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

def generateLANL(massR, chargeR, a, N_gamma=2.5e5):
    N_pot = 5.9e22 
    alpha = 1.0 / 137

#    a = 4  # num of layers

    m_e = 0.00051
    m_pi = 0.135
    m_eta = 0.548

    c_pi = 0.115
    c_eta = c_pi/30

    branch_pi = 0.98
    branch_eta = 0.39

    ageos = np.loadtxt("ageo_new.txt")
    ageo10mR = np.array(ageos[:,1]) / 2e7
    ageo35mR = np.array(ageos[:,2]) / 2e7
    ageo60mR = np.array(ageos[:,3]) / 2e7
    ageo100mR = np.array(ageos[:,4]) / 2e7

    ageo10m, _ = np.meshgrid(ageo10mR, chargeR)
    ageo35m, _ = np.meshgrid(ageo35mR, chargeR)
    ageo60m, _ = np.meshgrid(ageo60mR, chargeR)
    ageo100m, _ = np.meshgrid(ageo100mR, chargeR)

    mass, charge = np.meshgrid(massR, chargeR)

    def I_3(z, x):
        return 2 / (3 * 3.14) * ((1 - 4 * x / z) ** 0.5) * ((1 - z) ** 3) * (2 * x + z) / z ** 2

    def dl(x):
        return quad(I_3, 4 * x, 1, args=(x))[0]  # pass x to I
    v_dl = np.vectorize(dl)

    pi = np.zeros(np.shape(mass))
    eta = np.zeros(np.shape(mass))

    print(np.shape(mass))
    print(np.shape(ageo10m))
    print(mass)

    pi[mass < m_pi / 2] = ageo10m[mass < m_pi/2] * (1 - np.exp(-N_gamma * charge[mass < m_pi / 2] ** 2)) ** a * 2 * c_pi * branch_pi * alpha * N_pot * ageo10m[mass < m_pi/2] * v_dl(mass[mass < m_pi / 2] ** 2 / m_pi ** 2) * charge[mass < m_pi / 2] ** 2

    eta[mass < m_eta / 2] = ageo10m[mass < m_eta/2] * (1 - np.exp(-N_gamma * charge[mass < m_eta / 2] ** 2)) ** a * 2 * c_eta * branch_eta * alpha * N_pot * ageo10m[mass < m_eta/2] * v_dl(mass[mass < m_eta / 2] ** 2 / m_eta ** 2) * charge[mass < m_eta / 2] ** 2
    
    sensitivity10m = pi + eta

    pi = np.zeros(np.shape(mass))
    eta = np.zeros(np.shape(mass))

    pi[mass < m_pi / 2] = ageo35m[mass < m_pi/2] * (1 - np.exp(-N_gamma * charge[mass < m_pi / 2] ** 2)) ** a * 2 * c_pi * branch_pi * alpha * N_pot * ageo35m[mass < m_pi/2] * v_dl(mass[mass < m_pi / 2] ** 2 / m_pi ** 2) * charge[mass < m_pi / 2] ** 2

    eta[mass < m_eta / 2] = ageo35m[mass < m_eta/2] * (1 - np.exp(-N_gamma * charge[mass < m_eta / 2] ** 2)) ** a * 2 * c_eta * branch_eta * alpha * N_pot * ageo35m[mass < m_eta/2] * v_dl(mass[mass < m_eta / 2] ** 2 / m_eta ** 2) * charge[mass < m_eta / 2] ** 2

    sensitivity35m = pi + eta

    pi = np.zeros(np.shape(mass))
    eta = np.zeros(np.shape(mass))

    pi[mass < m_pi / 2] = ageo60m[mass < m_pi/2] * (1 - np.exp(-N_gamma * charge[mass < m_pi / 2] ** 2)) ** a * 2 * c_pi * branch_pi * alpha * N_pot * ageo35m[mass < m_pi/2] * v_dl(mass[mass < m_pi / 2] ** 2 / m_pi ** 2) * charge[mass < m_pi / 2] ** 2

    eta[mass < m_eta / 2] = ageo60m[mass < m_eta/2] * (1 - np.exp(-N_gamma * charge[mass < m_eta / 2] ** 2)) ** a * 2 * c_eta * branch_eta * alpha * N_pot * ageo35m[mass < m_eta/2] * v_dl(mass[mass < m_eta / 2] ** 2 / m_eta ** 2) * charge[mass < m_eta / 2] ** 2

    sensitivity60m = pi + eta

    pi = np.zeros(np.shape(mass))
    eta = np.zeros(np.shape(mass))

    pi[mass < m_pi / 2] = ageo100m[mass < m_pi/2] * (1 - np.exp(-N_gamma * charge[mass < m_pi / 2] ** 2)) ** a * 2 * c_pi * branch_pi * alpha * N_pot * ageo35m[mass < m_pi/2] * v_dl(mass[mass < m_pi / 2] ** 2 / m_pi ** 2) * charge[mass < m_pi / 2] ** 2

    eta[mass < m_eta / 2] = ageo100m[mass < m_eta/2] * (1 - np.exp(-N_gamma * charge[mass < m_eta / 2] ** 2)) ** a * 2 * c_eta * branch_eta * alpha * N_pot * ageo35m[mass < m_eta/2] * v_dl(mass[mass < m_eta / 2] ** 2 / m_eta ** 2) * charge[mass < m_eta / 2] ** 2

    sensitivity100m = pi + eta

    return sensitivity10m, sensitivity35m, sensitivity60m, sensitivity100m



if __name__ == '__main__':
    fig, ax = plt.subplots()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]')
    plt.ylabel('$\epsilon=Q/e$')
    plt.xlim(0.01, 10)
    plt.ylim(0.00001, 1)
    mass   = np.logspace(-2, 1, 500)
    charge = np.logspace(-5, 1, 500)
    masses, charges = np.meshgrid(mass, charge)
    print(np.shape(masses))
    submet = generateSubmet(masses, charges)
    submetctr = ax.contour(mass, charge, submet, levels=[20], colors = 'red')
    
    data = np.loadtxt("ageo_new.txt")
    mass = np.array(data[:,0]) / 1000
    lanl10m_2layer, lanl35m_2layer, lanl60m_2layer, lanl100m_2layer = generateLANL(mass, charge, 2)
    lanl10m_1layer, lanl35m_1layer, lanl60m_1layer, lanl100m_1layer = generateLANL(mass, charge, 1, 1.25e6)

    print("Generation completed") 



    lanl210ctr  = ax.contour(mass, charge, lanl10m_2layer, levels=[61], colors = 'skyblue')
    lanl235ctr  = ax.contour(mass, charge, lanl35m_2layer, levels=[61], colors = 'steelblue')
    lanl260ctr  = ax.contour(mass, charge, lanl60m_2layer, levels=[61], colors = 'deepskyblue')
    lanl2100ctr  = ax.contour(mass, charge, lanl100m_2layer, levels=[61], colors = 'lightskyblue')
    lanl110ctr = ax.contour(mass, charge, lanl10m_1layer, levels=[428], colors='blueviolet')
    lanl135ctr = ax.contour(mass, charge, lanl35m_1layer, levels=[428], colors='darkviolet')
    lanl160ctr = ax.contour(mass, charge, lanl60m_1layer, levels=[428], colors='violet')
    lanl1100ctr = ax.contour(mass, charge, lanl100m_1layer, levels=[428], colors='purple')

    h1, _ = submetctr.legend_elements()
    h2, _ = lanl210ctr.legend_elements()
    h3, _ = lanl235ctr.legend_elements()
    h4, _ = lanl260ctr.legend_elements()
    h5, _ = lanl2100ctr.legend_elements()
    h6, _ = lanl110ctr.legend_elements()
    h7, _ = lanl135ctr.legend_elements()
    h8, _ = lanl160ctr.legend_elements()
    h9, _ = lanl1100ctr.legend_elements()


    ax.legend([h1[0], h2[0], h3[0], h4[0], h5[0], h6[0], h7[0], h9[0]], 
            ['SUBMET', 
            '2 Layers, 80 bars, bkg = 931, 10m', '2 Layers, 80 bars,  bkg = 931, 35m', '2 Layers, 80 bars,  bkg = 931, 60m', '2 Layers, 80 bars,  bkg = 931, 100m', 
            '1 Layers, 1 bar  bkg = 47k, 5Npe, 10m', '1 Layers, 1 bar,  bkg = 47k, 5Npe, 35m', '1 Layer, 1 bar, bkg = 47k, 5Npe, 60m', '1 Layer, 1 bar, bkg = 47k, 5Npe, 100m'], loc = 'lower right')

    plt.savefig('sensitivity_lanl_new.pdf')
    plt.show()
