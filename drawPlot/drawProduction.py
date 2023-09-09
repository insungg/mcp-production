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


def mesonDecayProduction(fileName, NPOT):
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

    # read geometric acceptance obtained from PYTHIA
#    mass = np.array([])
#    ageo_pi = np.array([])
#    ageo_eta = np.array([])
#    ageo_jpsi = np.array([])
#    ageo_upsilon = np.array([])
# with open(fileName, 'r') as f:
#         data = f.readlines()
#         for line in data:
#             values = line.strip().split()
#             np.append(mass,         float(values[0]))
#             np.append(ageo_pi,      float(values[1]))
#             np.append(ageo_eta,     float(values[2]))
#             np.append(ageo_jpsi,    float(values[3]))
#             np.append(ageo_upsilon, float(values[3])) # use same ageo with jpsi
    ageo_data = np.loadtxt(fileName)
    mass        = ageo_data[:, 0]
    ageo_pi     = ageo_data[:, 1]
    ageo_eta    = ageo_data[:, 2]
    ageo_jpsi   = ageo_data[:, 3]
    ageo_upsilon = ageo_jpsi
                
# ageo_jpsi.fill(1e-3)
    ageo_upsilon.fill(1e-3)
    print(ageo_jpsi)
    print(ageo_upsilon)

    # define phase space integrals 
    def I2(x, y):
        return ((1 + 2 * x) * np.sqrt(1 - 4 * x)) / ((1 + 2 * y) * np.sqrt(1 - 4 * y))

    def I3_integrand(z, x):
        return 2 / (3 * 3.14) * np.sqrt(1 - 4 * x / z)  * ((1 - z) ** 3) * (2 * x + z) / (z ** 2)

    def I3(x):
        return quad(I3_integrand, 4 * x, 1, args=(x))[0]  # integrate
    v_I3 = np.vectorize(I3) # vectorize

    # define productions
    pi      = np.zeros(np.shape(mass))
    eta     = np.zeros(np.shape(mass))
    jpsi    = np.zeros(np.shape(mass))
    upsilon = np.zeros(np.shape(mass))

    # masking to achieve kinematical validity
    # Dalitz decay
#    pi[mass] = NPOT * ageo_pi[mass] * 2 * c_pi * branch_pi * alpha * v_I3( mass[mass] ** 2 / mass[mass] ** 2)
#    eta[mass] = NPOT * ageo_eta[mass] * 2 * c_eta * branch_eta * alpha * v_I3( mass[mass] ** 2 / mass[mass] ** 2)
#    # Direct decay
#    jpsi[mass] = NPOT * ageo_jpsi[mass] * 2  * c_jpsi * branch_jpsi * I2( mass[mass] ** 2 / m_jpsi ** 2, m_e ** 2 / m_jpsi ** 2  )
#    upsilon[mass] = NPOT * ageo_upsilon[mass] * 2 * c_upsilon * branch_upsilon * I2( mass[mass] ** 2 / m_upsilon ** 2, m_e ** 2 / m_upsilon ** 2  )
 
    # Dalitz decay
    pi[mass < m_pi/2] = NPOT * ageo_pi[mass < m_pi/2] * 2 * c_pi * branch_pi * alpha * v_I3( mass[mass < m_pi/2] ** 2 / m_pi ** 2)
    eta[mass < m_eta/2] = NPOT * ageo_eta[mass < m_eta/2] * 2 * c_eta * branch_eta * alpha * v_I3( mass[mass < m_eta/2] ** 2 /  m_eta/2 ** 2)
    # Direct decay
    jpsi[mass < m_jpsi/2] = NPOT * ageo_jpsi[mass < m_jpsi/2] * 2  * c_jpsi * branch_jpsi * I2( mass[mass < m_jpsi/2] ** 2 / m_jpsi ** 2, m_e ** 2 / m_jpsi ** 2  )
    upsilon[mass < m_upsilon/2] = NPOT * ageo_upsilon[mass < m_upsilon/2] * 2 * c_upsilon * branch_upsilon * I2( mass[mass < m_upsilon/2] ** 2 / m_upsilon ** 2, m_e ** 2 / m_upsilon ** 2  )

    return pi, eta, jpsi, upsilon


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
    mass     = ageo_data[:, 0]
    cross_dy = ageo_data[:, 1] * 1e-12 # pico barn
    ageo_dy  = ageo_data[:, 2]

    return NPOT * cross_dy / totalCrossSection * ageo_dy



if __name__ == '__main__':
    # read the original mass file
#    mass = np.array([])
#    with open('mass.txt', 'r') as f:
#        data = f.readline()
#        for line in data:
# values = line.strip().split()
#            print(line)
#            np.append(mass, float(line) )
    mass = np.loadtxt('mass.txt')

    # get decay production
    numi = 6e20
    dune = 1e21
#   pi_numi, eta_numi, jpsi_numi, upsilon_numi = mesonDecayProduction(numi)
#   pi_dune, eta_dune, jpsi_dune, upsilon_dune = mesonDecayProduction(dune)

    mesonDecay_numi = list(mesonDecayProduction('ageo_decay_numi.txt', numi))
    mesonDecay_dune = list(mesonDecayProduction('ageo_decay_dune.txt', dune))

    # get DY production
    dy_numi = dyProduction('ageo_dy_numi.txt', numi)
    dy_dune = dyProduction('ageo_dy_dune.txt', dune)

    # add padings
    for i in range(len(mesonDecay_numi)):
        difference = mass.size - mesonDecay_numi[i].size
        mesonDecay_numi[i] = np.pad(mesonDecay_numi[i], (0, difference), 'constant',  constant_values=0)
    for i in range(len(mesonDecay_dune)):
        difference = mass.size - mesonDecay_dune[i].size
        mesonDecay_dune[i] = np.pad(mesonDecay_dune[i], (0, difference), 'constant', constant_values=0)

    difference = mass.size - dy_numi.size
    dy_numi = np.pad(dy_numi, (0, difference), 'constant', constant_values=0)
    difference = mass.size - dy_dune.size 
    dy_dune = np.pad(dy_dune, (0, difference), 'constant', constant_values=0)

    total_numi = dy_numi
    for i in range(len(mesonDecay_numi)):
        total_numi = total_numi + mesonDecay_numi[i]
    
    total_dune = dy_dune
    for i in range(len(mesonDecay_dune)):
        total_dune = total_dune + mesonDecay_dune[i]
    print(mass.size)
    print(dy_numi.size)
    print(mesonDecay_numi[1].size)
    

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
    ax1.plot(mass, mesonDecay_numi[0], label = r'$\pi^{0}\to\gamma\chi\overline{\chi}$')
    ax1.plot(mass, mesonDecay_numi[1], label = r'$\eta\to\gamma\chi\overline{\chi}$')
    ax1.plot(mass, mesonDecay_numi[2], label = r'$J/\psi\to\chi\overline{\chi}$')
    ax1.plot(mass, mesonDecay_numi[3], label = r'$\Upsilon\to\chi\overline{\chi}$')
    ax1.plot(mass, dy_numi,            label =  r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$')
    ax1.plot(mass, total_numi,         label = 'Total')
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
    ax2.plot(mass, mesonDecay_dune[0], label = r'$\pi^{0}\to\gamma\chi\overline{\chi}$')
    ax2.plot(mass, mesonDecay_dune[1], label = r'$\eta\to\gamma\chi\overline{\chi}$')
    ax2.plot(mass, mesonDecay_dune[2], label = r'$J/\psi\to\chi\overline{\chi}$')
    ax2.plot(mass, mesonDecay_dune[3], label = r'$\Upsilon\to\chi\overline{\chi}$')
    ax2.plot(mass, dy_dune,            label = r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$')
    ax2.plot(mass, total_dune,         label = 'Total')
    ax2.legend(loc = 'lower left', ncol = 2)
    ax2.text(0.55, 0.93, '$1$ Year at DUNE ($10^{21}$ POT) \n $574$ m, $1$ m $\\times$ $1$ m detector', transform=ax2.transAxes, fontsize=10, verticalalignment='center', multialignment = 'left', bbox=textbox_props)

    plt.show()

    # save plot as pgf format
    fig.savefig('MCP-production.pgf')

