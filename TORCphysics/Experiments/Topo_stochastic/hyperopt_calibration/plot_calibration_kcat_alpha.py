from TORCphysics import Circuit
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
from TORCphysics import visualization as vs
import matplotlib.pyplot as plt
import numpy as np

# Initialise figure and circuit initial conditions
# ----------------------------------------------------------------------------------------------------------------------

# Some initial conditions
nsim = 100  # number of simulations per test
topo_k_on_0 = .005  # 0.0075
topo_k_cat_0 = 10.0
topo_concentration_0 = .25
gyra_k_on_0 = .005  # 0.0075
gyra_k_cat_0 = -20
gyra_concentration_0 = 0.0

# Figure
width = 8
height = 4
figure_output = 'calibration_kcats-alpha.png'

fig, axs = plt.subplots(4, figsize=(width, 4 * height), tight_layout=True)

# Circuit conditions
circuit_filename_topoI_0 = 'circuit_3000bp_negative.csv'
circuit_filename_gyrase_0 = 'circuit_3000bp_positive.csv'
sites_filename_0 = 'sites.csv'
enzymes_filename_0 = 'enzymes.csv'
environment_continuum_filename_0 = 'environment_continuum.csv'
environment_stochastic_filename_0 = 'environment_stochastic.csv'
output_prefix = 'output'
frames = 2000
series = True
continuation = False
dt = .5

ctopo_k_cat = 'topoI_kcat-alpha.txt'
my_topo = np.loadtxt(ctopo_k_cat)
my_topo_k_cat =my_topo[0]
my_topo_alpha =my_topo[1]
cgyra_k_cat = 'gyra_kcat-alpha.txt'
my_gyra = np.loadtxt(cgyra_k_cat)
my_gyra_k_cat =my_gyra[0]
my_gyra_alpha =my_gyra[1]


# Some functions
# ----------------------------------------------------------------------------------------------------------------------
def run_many_stochastic(topo_concentration, topo_k_on, topo_k_cat, gyra_concentration, gyra_k_on, gyra_k_cat,
                        circuit_filename, sites_filename, enzymes_filename, environment_filename,
                        output_prefix, nframes, series, continuation, dt, tm, mm):
    supercoiling = np.zeros((nframes + 1, nsim))
    for i in range(nsim):
        supercoiling[:, i] = run_stochastic_sim(topo_concentration, topo_k_on, topo_k_cat,
                                                gyra_concentration, gyra_k_on, gyra_k_cat,
                                                circuit_filename, sites_filename, enzymes_filename,
                                                environment_filename, output_prefix, nframes, series, continuation,
                                                dt, tm, mm)
    return supercoiling


def run_stochastic_sim(topo_concentration, topo_k_on, topo_k_cat, gyra_concentration, gyra_k_on, gyra_k_cat,
                       circuit_filename, sites_filename, enzymes_filename, environment_filename,
                       output_prefix, nframes, series, continuation, dt, tm, mm):
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, nframes, series, continuation, dt, tm, mm)
    my_circuit.environmental_list[0].k_on = topo_k_on
    my_circuit.environmental_list[0].k_cat = topo_k_cat
    my_circuit.environmental_list[0].concentration = topo_concentration
    my_circuit.environmental_list[1].k_on = gyra_k_on
    my_circuit.environmental_list[1].k_cat = gyra_k_cat
    my_circuit.environmental_list[1].concentration = gyra_concentration

    my_supercoiling = np.zeros(nframes + 1)
    my_supercoiling[0] = my_circuit.superhelical

    # run simulation
    for frame in range(1, nframes + 1):
        # BINDING
        # --------------------------------------------------------------
        new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                           my_circuit.rng)
        my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

        # EFFECT
        # --------------------------------------------------------------
        effects_list = em.effect_model(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                       my_circuit.topoisomerase_model, my_circuit.mechanical_model)
        my_circuit.apply_effects(effects_list)

        # UNBINDING
        # --------------------------------------------------------------
        drop_list_index, drop_list_enzyme = bm.unbinding_model(my_circuit.enzyme_list, my_circuit.dt,
                                                               my_circuit.rng)
        my_circuit.drop_enzymes(drop_list_index)
        my_circuit.add_to_environment(drop_list_enzyme)
        # UPDATE GLOBALS
        # --------------------------------------------------------------
        my_circuit.update_global_twist()
        my_circuit.update_global_superhelical()
        my_supercoiling[frame] = my_circuit.superhelical
    return my_supercoiling


# Plot recognition curves
# ----------------------------------------------------------------------------------------------------------------------
ax = axs[3]
circuit_filename_0 = circuit_filename_topoI_0
# Load circuits
continuum_circuit = Circuit(circuit_filename_0, sites_filename_0, enzymes_filename_0, environment_continuum_filename_0,
                            output_prefix, frames, series, continuation, dt, 'continuum', 'uniform')
stochastic_circuit = Circuit(circuit_filename_0, sites_filename_0, enzymes_filename_0,
                             environment_stochastic_filename_0,
                             output_prefix, frames, series, continuation, dt, 'stochastic', 'uniform')
stochastic_circuit.environmental_list[0].k_cat = my_topo_k_cat
stochastic_circuit.environmental_list[1].k_cat = my_gyra_k_cat
stochastic_circuit.environmental_list[0].concentration = \
    my_topo_alpha * stochastic_circuit.environmental_list[0].concentration
stochastic_circuit.environmental_list[1].concentration = \
    my_gyra_alpha * stochastic_circuit.environmental_list[1].concentration
continuum_circuit.environmental_list[0].k_cat = my_topo_k_cat
continuum_circuit.environmental_list[1].k_cat = my_gyra_k_cat
topo_effective_concentration = stochastic_circuit.environmental_list[0].concentration
gyra_effective_concentration = stochastic_circuit.environmental_list[1].concentration

vs.plot_site_response_curves(stochastic_circuit, ax)
ax.set_title('Topoisomerases response curves')

# Plot global supercoiling responses - topo_I calibration
# ---------------------------------------------------------
ax = axs[0]
circuit_filename_0 = circuit_filename_topoI_0
topo_concentration_0 = 0.25
gyra_concentration_0 = 0.0
time = np.arange(0, (frames + 1) * dt, dt)

continuum_circuit = Circuit(circuit_filename_0, sites_filename_0, enzymes_filename_0, environment_continuum_filename_0,
                            output_prefix, frames, series, continuation, dt, 'continuum', 'uniform')
continuum_circuit.environmental_list[1].concentration = gyra_concentration_0  # Turn gyrase concentration to 0
# Run continuum case and get its global sigma
continuum_circuit.run()
mask = continuum_circuit.sites_df['type'] == 'circuit'
sigma_continuum = continuum_circuit.sites_df[mask]['superhelical'].to_numpy()

# And run the stochastic case
sigma_topoI = run_many_stochastic(topo_effective_concentration, topo_k_on_0, my_topo_k_cat,
                                        gyra_concentration_0, gyra_k_on_0, gyra_k_cat_0,
                                        circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                        environment_stochastic_filename_0, output_prefix, frames, series, continuation,
                                        dt, 'stochastic', 'uniform')
sigma_topoI = np.mean(sigma_topoI, axis=1)

ax.plot(time, sigma_continuum, 'black', label='continuum')
ax.plot(time, sigma_topoI, 'red', label='stochastic')
ax.grid(True)
ax.set_xlabel(r'time ($s$)')
ax.set_ylabel(r'$\sigma$')
ax.set_title('Topoisomerase I calibration')
ax.legend(loc='best')

#k_cat
text1 = r'$k_{cat}=$'
text2 = f'{my_topo_k_cat:.2f}'
text3 = r'$s^{-1}$'
text = text1 + text2 + text3
xs = 0.78 #500
ys = 0.55 # -0.02
props = dict(boxstyle='round', facecolor='gray', alpha=0.4)
ax.text(xs, ys, text, transform=ax.transAxes, verticalalignment='top', bbox=props)

#alpha
text1 = r'$\alpha=$'
text2 = f'{my_topo_alpha:.2f}'
text3 = r'$\mu M^{-1}$'
text = text1 + text2 + text3
xs = 0.78 #500
ys = 0.45 # -0.02
props = dict(boxstyle='round', facecolor='gray', alpha=0.4)
ax.text(xs, ys, text, transform=ax.transAxes, verticalalignment='top', bbox=props)

# Plot global supercoiling responses - gyrase calibration
# ---------------------------------------------------------
ax = axs[1]
circuit_filename_0 = circuit_filename_gyrase_0
topo_concentration_0 = 0.0
gyra_concentration_0 = 0.25
time = np.arange(0, (frames + 1) * dt, dt)

continuum_circuit = Circuit(circuit_filename_0, sites_filename_0, enzymes_filename_0, environment_continuum_filename_0,
                            output_prefix, frames, series, continuation, dt, 'continuum', 'uniform')
continuum_circuit.environmental_list[0].concentration = topo_concentration_0  # Turn topo I concentration to 0
# Run continuum case and get its global sigma
continuum_circuit.run()
mask = continuum_circuit.sites_df['type'] == 'circuit'
sigma_continuum = continuum_circuit.sites_df[mask]['superhelical'].to_numpy()

# And run the stochastic case
sigma_topoI = run_many_stochastic(topo_concentration_0, topo_k_on_0, my_topo_k_cat,
                                        gyra_effective_concentration, gyra_k_on_0, my_gyra_k_cat,
                                        circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                        environment_stochastic_filename_0, output_prefix, frames, series, continuation,
                                        dt, 'stochastic', 'uniform')
sigma_topoI = np.mean(sigma_topoI, axis=1)

ax.plot(time, sigma_continuum, 'black', label='continuum')
ax.plot(time, sigma_topoI, 'red', label='stochastic')
ax.grid(True)
ax.set_xlabel(r'time ($s$)')
ax.set_ylabel(r'$\sigma$')
ax.set_title('Gyrase calibration')
ax.legend(loc='best')

#k_cat
text1 = r'$k_{cat}=$'
text2 = f'{my_gyra_k_cat:.2f}'
text3 = r'$s^{-1}$'
text = text1 + text2 + text3
xs = 0.6 #500
ys = 0.95  # -0.02
props = dict(boxstyle='round', facecolor='gray', alpha=0.4)
ax.text(xs, ys, text, transform=ax.transAxes, verticalalignment='top', bbox=props)

#alpha
text1 = r'$\alpha=$'
text2 = f'{my_gyra_alpha:.2f}'
text3 = r'$\mu M^{-1}$'
text = text1 + text2 + text3
xs = 0.6 #500
ys = 0.85 # -0.02
props = dict(boxstyle='round', facecolor='gray', alpha=0.4)
ax.text(xs, ys, text, transform=ax.transAxes, verticalalignment='top', bbox=props)

# Plot global supercoiling responses - both enzymes active
# ---------------------------------------------------------
ax = axs[2]
circuit_filename_0 = circuit_filename_topoI_0
topo_concentration_0 = 0.25
gyra_concentration_0 = 0.25
time = np.arange(0, (frames + 1) * dt, dt)

continuum_circuit = Circuit(circuit_filename_0, sites_filename_0, enzymes_filename_0, environment_continuum_filename_0,
                            output_prefix, frames, series, continuation, dt, 'continuum', 'uniform')
# Run continuum case and get its global sigma
continuum_circuit.run()
mask = continuum_circuit.sites_df['type'] == 'circuit'
sigma_continuum = continuum_circuit.sites_df[mask]['superhelical'].to_numpy()

# And run the stochastic case
sigma_topoI = run_many_stochastic(topo_effective_concentration, topo_k_on_0, my_topo_k_cat,
                                        gyra_effective_concentration, gyra_k_on_0, my_gyra_k_cat,
                                        circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                        environment_stochastic_filename_0, output_prefix, frames, series, continuation,
                                        dt, 'stochastic', 'uniform')
sigma_topoI = np.mean(sigma_topoI, axis=1)

ax.plot(time, sigma_continuum, 'black', label='continuum')
ax.plot(time, sigma_topoI, 'red', label='stochastic')
ax.grid(True)
ax.set_xlabel(r'time ($s$)')
ax.set_ylabel(r'$\sigma$')
ax.set_title('both enzymes active')
ax.legend(loc='best')

plt.savefig(figure_output)
