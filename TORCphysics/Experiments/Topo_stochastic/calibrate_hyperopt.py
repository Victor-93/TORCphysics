from TORCphysics import Circuit
import matplotlib.pyplot as plt
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import sys
from hyperopt import tpe, hp, fmin
import numpy as np

# Now we will calibrate by minimizing an objective function using hyperopt


# System conditions
lengths = [3000, 4000]  # , 5000, 6000]
nsim = 10

# Some initial conditions
topo_k_on_0 = .005  # 0.0075
topo_k_cat_0 = 10.0
topo_concentration_0 = .25
gyra_k_on_0 = .005  # 0.0075
gyra_k_cat_0 = -20
gyra_concentration_0 = 0.0
superhelical_topo = -0.04  # Initial condition for calibrating topo
superhelical_gyrase = 0.04  # And for gyrase

# Circuit conditions
circuit_filename_0 = 'circuit.csv'
sites_filename_0 = 'sites.csv'
enzymes_filename_0 = 'enzymes.csv'
environment_continuum_filename_0 = 'environment_continuum.csv'
environment_filename_0 = 'environment.csv'
tm = 'stochastic'
output_prefix = 'continuum'
frames = 1000
series = True
continuation = False
mm = 'uniform'
dt = .5

time = np.arange(0, (frames + 1) * dt, dt)


def objective_topoI_k_cat(params):
    my_topo_k_cat = params['topo_k_cat']
    supercoiling = run_many_stochastic(topo_concentration_0, topo_k_on_0, my_topo_k_cat,
                                       gyra_concentration_0, gyra_k_on_0, gyra_k_cat_0,
                                       circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                       environment_filename_0, output_prefix, frames, series, continuation,
                                       dt, 'stochastic', mm)
    my_objective = np.sum(np.square(np.mean(supercoiling, axis=1) - sigma_continuum))
    return my_objective


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


# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# I will be plotting as the experiment advances.

# Initialise figure and circuit initial conditions
# ----------------------------------------------------------------------------------------------------------------------

# Figure
width = 8
height = 3
figure_output = 'calibration_hyperopt.png'

fig, axs = plt.subplots(1, figsize=(width, 1.25 * height), tight_layout=True)

# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------
tm = 'continuum'
continuum_circuit = Circuit(circuit_filename_0, sites_filename_0, enzymes_filename_0, environment_continuum_filename_0,
                            output_prefix, frames, series, continuation, dt, tm, mm)
continuum_circuit.environmental_list[1].concentration = gyra_concentration_0
continuum_circuit.run()

# Get global supercoiling responses from the continuum case
mask = continuum_circuit.sites_df['type'] == 'circuit'  # This one contains global superhelical density
sigma_continuum = continuum_circuit.sites_df[mask]['superhelical'].to_numpy()

ax = axs  # [0]
ax.plot(time, sigma_continuum, color='black')

# Stochastic 0
# ----------------------------------------------------------------------------------------------------------------------
supercoiling_0 = run_many_stochastic(topo_concentration_0, topo_k_on_0, topo_k_cat_0,
                                     gyra_concentration_0, gyra_k_on_0, gyra_k_cat_0,
                                     circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                     environment_filename_0, output_prefix, frames, series, continuation,
                                     dt, 'stochastic', mm)

supercoiling_0 = np.mean(supercoiling_0, axis=1)
ax.plot(time, supercoiling_0, color='green')

# Optimization case
# ----------------------------------------------------------------------------------------------------------------------
space = {
    'topo_k_cat': hp.uniform('topo_k_cat', 5.0, 30.0)
}
best = fmin(
    fn=objective_topoI_k_cat, # Objective Function to optimize
    space=space,  # Hyperparameter's Search Space
    algo=tpe.suggest, # Optimization algorithm (representative TPE)
    max_evals=100 # Number of optimization attempts
)
print(best)
print(best.values())
topo_k_cat_best = best['topo_k_cat']
supercoiling_best = run_many_stochastic(topo_concentration_0, topo_k_on_0, topo_k_cat_best,
                                     gyra_concentration_0, gyra_k_on_0, gyra_k_cat_0,
                                     circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                     environment_filename_0, output_prefix, frames, series, continuation,
                                     dt, 'stochastic', mm)

supercoiling_best = np.mean(supercoiling_best, axis=1)
ax.plot(time, supercoiling_best, color='red')

plt.savefig(figure_output)

