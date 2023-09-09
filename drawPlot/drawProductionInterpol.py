import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter

# define matplotlib pgf output settings
# matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})


from scipy.integrate import quad

# define matplotlib parameters
# SMALL_SIZE = 20
# MEDIUM_SIZE = 20
# BIGGER_SIZE = 20
# 
# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE + 5)  # fontsize of the figure title


def mesonDecayProductionLight(fileName, NPOT):
    # EM constant
    alpha = 1.0 / 137

    # mass in GeV 
    m_e = 0.00051
    m_pi = 0.135
    m_eta = 0.548
    # m_rho = 0.775
    # m_omega = 0.782
    # m_phi = 1.019
    m_jpsi = 3.096
    m_upsilon = 9.46

    # meson / NPOT obtained from PYTHIA
    c_pi = 3.98
    c_eta = 0.51
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

    ageo_data   = np.loadtxt(fileName)
    mass        = np.array(ageo_data[:, 0])
    ageo_pi     = np.array(ageo_data[:, 1])
    ageo_eta    = np.array(ageo_data[:, 2])
                
    # define phase space integrals 
    def I3_integrand(z, x):
        return 2 / (3 * 3.14) * np.sqrt(1 - 4 * x / z)  * ((1 - z) ** 3) * (2 * x + z) / (z ** 2)

    def I3(x):
        return quad(I3_integrand, 4 * x, 1, args=(x))[0]  # integrate
    v_I3 = np.vectorize(I3) # vectorize

    # define productions
    pi      = np.zeros(np.shape(mass))
    eta     = np.zeros(np.shape(mass))
 
    # Dalitz decay
    pi[mass < m_pi/2] = NPOT * ageo_pi[mass < m_pi/2] * 2 * c_pi * branch_pi * alpha * v_I3( mass[mass < m_pi/2] ** 2 / m_pi ** 2)
    eta[mass < m_eta/2] = NPOT * ageo_eta[mass < m_eta/2] * 2 * c_eta * branch_eta * alpha * v_I3( mass[mass < m_eta/2] ** 2 /  m_eta ** 2)

    return mass, pi, eta

def mesonDecayProductionHeavy(fileName, NPOT):
    # EM constant
    alpha = 1.0 / 137

    # mass in GeV 
    m_e = 0.00051
    m_pi = 0.135
    m_eta = 0.548
    # m_rho = 0.775
    # m_omega = 0.782
    # m_phi = 1.019
    m_jpsi = 3.096
    m_upsilon = 9.46

    # meson / NPOT obtained from PYTHIA
    c_pi = 3.98
    c_eta = 0.51
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

    ageo_data   = np.loadtxt(fileName)
    mass        = np.array(ageo_data[:, 0])
    ageo_jpsi   = np.array(ageo_data[:, 1])
    ageo_upsilon = ageo_jpsi

    ageo_upsilon.fill(1e-3)
                
    # define phase space integrals 
    def I2(x, y):
        return ((1 + 2 * x) * np.sqrt(1 - 4 * x)) / ((1 + 2 * y) * np.sqrt(1 - 4 * y))


    # define productions
    jpsi    = np.zeros(np.shape(mass))
    upsilon = np.zeros(np.shape(mass))

    # Direct decay
    jpsi[mass < m_jpsi/2] = NPOT * ageo_jpsi[mass < m_jpsi/2] * 2  * c_jpsi * branch_jpsi * I2( mass[mass < m_jpsi/2] ** 2 / m_jpsi ** 2, m_e ** 2 / m_jpsi ** 2  )
    upsilon[mass < m_upsilon/2] = NPOT * ageo_upsilon[mass < m_upsilon/2] * 2 * c_upsilon * branch_upsilon * I2( mass[mass < m_upsilon/2] ** 2 / m_upsilon ** 2, m_e ** 2 / m_upsilon ** 2  )

    return mass, jpsi, upsilon

def dyProduction(fileName, NPOT):
    totalCrossSection = 300e-3

#     mass = np.array([])
#     cross_dy = np.array([])
#     ageo_dy = np.array([])
# 
#     with open(fileName, 'r') as f:
#         data = f.readlines()
#         for line in data:
#             values = line.strip().split()
#             np.append(cross_dy, float(values[1]) * 1e-12) # convert to pico-barn
#             np.append(ageo_dy,  float(values[2]))

    ageo_data = np.loadtxt(fileName)
    mass     = np.array(ageo_data[:, 0])
    cross_dy = np.array(ageo_data[:, 1] * 1e-12) # pico barn
    ageo_dy  = np.array(ageo_data[:, 2])

    return mass, NPOT * cross_dy / totalCrossSection * ageo_dy



if __name__ == '__main__':
    # get decay production
    numi = 6e20
    dune = 1e21

    # read data and generate production
    mass_numi_light, pi_numi, eta_numi       = list(mesonDecayProductionLight('ageo_decay_numi_light.txt', numi))
    mass_numi_heavy, jpsi_numi, upsilon_numi = list(mesonDecayProductionHeavy('ageo_decay_numi_heavy.txt', numi))
    mass_dune_light, pi_dune, eta_dune       = list(mesonDecayProductionLight('ageo_decay_dune_light.txt', dune))
    mass_dune_heavy, jpsi_dune, upsilon_dune = list(mesonDecayProductionHeavy('ageo_decay_dune_heavy.txt', dune))

    # get DY production
    mass_numi_dy, dy_numi = dyProduction('ageo_dy_numi_full.txt', numi)
    mass_dune_dy, dy_dune = dyProduction('ageo_dy_dune_full.txt', dune)

    # get a sorted list of all mass values
    mass_numi_total = np.unique(np.concatenate((np.concatenate((mass_numi_light, mass_numi_heavy), axis=None), mass_numi_dy), axis=None))
    mass_dune_total = np.unique(np.concatenate((np.concatenate((mass_dune_light, mass_dune_heavy), axis=None), mass_dune_dy), axis=None))

    # interpolate each contributions to obtain total acceptance
    pi_numi_i       = np.interp(mass_numi_total, mass_numi_light, pi_numi,      left=0, right=0)
    eta_numi_i      = np.interp(mass_numi_total, mass_numi_light, eta_numi,     left=0, right=0)
    jpsi_numi_i     = np.interp(mass_numi_total, mass_numi_heavy, jpsi_numi,    left=0, right=0)
    upsilon_numi_i  = np.interp(mass_numi_total, mass_numi_heavy, upsilon_numi, left=0, right=0)
    dy_numi_i       = np.interp(mass_numi_total, mass_numi_dy,    dy_numi,      left=0, right=0)

    total_numi = pi_numi_i + eta_numi_i + jpsi_numi_i + upsilon_numi_i + dy_numi_i

    pi_dune_i       = np.interp(mass_dune_total, mass_dune_light, pi_dune,      left=0, right=0)
    eta_dune_i      = np.interp(mass_dune_total, mass_dune_light, eta_dune,     left=0, right=0)
    jpsi_dune_i     = np.interp(mass_dune_total, mass_dune_heavy, jpsi_dune,    left=0, right=0)
    upsilon_dune_i  = np.interp(mass_dune_total, mass_dune_heavy, upsilon_dune, left=0, right=0)
    dy_dune_i       = np.interp(mass_dune_total, mass_dune_dy,    dy_dune,      left=0, right=0)

    total_dune = pi_dune_i + eta_dune_i + jpsi_dune_i + upsilon_dune_i + dy_dune_i
    

    # draw plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    plt.subplots_adjust(left=0.08, right=0.92, wspace=0.3)
    x_ticks = [10**i for i in range(-2, 2)]
    x_tick_labels = ['$10^{-2}$', '$10^{-1}$', '1', '$10^{1}$']
    y_ticks = [10**i for i in range(0, 19)]
    y_tick_labels = ['$1$'] + [(lambda x : f'$10^{{{x}}}$' if x % 2 == 0 else ' ')(i) for i in range(1, 19)]
    textbox_props = dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.7)

    ax1.margins(x=0, y=0)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$m_{\chi}$ [GeV]')
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_tick_labels)
    ax1.set_ylabel(r'$N_{\chi} / \epsilon^{2}$')
    ax1.set_ylim(1, 1e18)
    ax1.set_yticks(y_ticks)
    ax1.set_yticklabels(y_tick_labels)
    ax1.plot(mass_numi_light, pi_numi, label = r'$\pi^{0}\to\gamma\chi\overline{\chi}$')
    ax1.plot(mass_numi_light, eta_numi, label = r'$\eta\to\gamma\chi\overline{\chi}$')
    ax1.plot(mass_numi_heavy, jpsi_numi, label = r'$J/\psi\to\chi\overline{\chi}$')
    ax1.plot(mass_numi_heavy, upsilon_numi, label = r'$\Upsilon\to\chi\overline{\chi}$')
    ax1.plot(mass_numi_dy, dy_numi, label =  r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$')
    ax1.plot(mass_numi_total, total_numi, label = 'Total')
    ax1.legend(loc = 'lower left', ncol = 2)
    ax1.text(0.55, 0.93, '$1$ Year at NuMI ($6\\times 10^{20}$ POT) \n $1040$ m, $1$ m $\\times$ $1$ m  detector', transform=ax1.transAxes, fontsize=10, verticalalignment='center', multialignment = 'left',  bbox=textbox_props)

    ax2.margins(x=0, y=0)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel(r'$m_{\chi}$ [GeV]')
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_tick_labels)
    ax2.set_ylabel(r'$N_{\chi} / \epsilon^{2}$')
    ax2.set_ylim(1, 1e18)
    ax2.set_yticks(y_ticks)
    ax2.set_yticklabels(y_tick_labels)
    ax2.plot(mass_dune_light, pi_dune, label = r'$\pi^{0}\to\gamma\chi\overline{\chi}$')
    ax2.plot(mass_dune_light, eta_dune, label = r'$\eta\to\gamma\chi\overline{\chi}$')
    ax2.plot(mass_dune_heavy, jpsi_dune, label = r'$J/\psi\to\chi\overline{\chi}$')
    ax2.plot(mass_dune_heavy, upsilon_dune, label = r'$\Upsilon\to\chi\overline{\chi}$')
    ax2.plot(mass_dune_dy, dy_dune, label = r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$')
    ax2.plot(mass_dune_total, total_dune,         label = 'Total')
    ax2.legend(loc = 'lower left', ncol = 2)
    ax2.text(0.55, 0.93, '$1$ Year at DUNE ($10^{21}$ POT) \n $574$ m, $1$ m $\\times$ $1$ m detector', transform=ax2.transAxes, fontsize=10, verticalalignment='center', multialignment = 'left', bbox=textbox_props)

    plt.show()

    # save plot as pgf format
    fig.savefig('MCP-production.pdf')

