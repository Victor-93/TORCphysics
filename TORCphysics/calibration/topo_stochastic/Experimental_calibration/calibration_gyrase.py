import sys

import matplotlib.pyplot as plt
import numpy as np
from hyperopt import tpe, hp, fmin
import multiprocessing
from TORCphysics import parallelization_tools as pt
# import pickle
import pandas as pd

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# TODO: Add a proper description
# TODO: Set the number of tests/simulations
# TODO: Save plot!
# This is just a test to reproduce the global supercoiling response curves from the paper:
# Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching

# We now have estimations of the initial and final superhelical densities.
# For now, let's not consider ATP.
# Let's run examples of the actual calibration. Remember that in your simulation, you need to consider the
# density of the plasmid concentration as well, so you have many concentrations to test.

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 0.5
initial_time = 0
final_time = 600
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# For the simulation
circuit_filename = 'circuit_test.csv'
sites_filename = 'sites_test.csv'
enzymes_filename = 'enzymes_test.csv'
environment_filename = 'environment_test.csv'
tm = 'stochastic'
output_prefix = 'test0'
series = True
continuation = False
mm = 'uniform'

# For parallelization and calibration
n_simulations = 5
tests = 2  # 00  # number of tests for parametrization

my_vars = ['k_cat', 'k_on', 'k_off', 'width', 'threshold']

# RANGES FOR RANDOM SEARCH
# -----------------------------------
# Gyrase ranges
gyrase_out = 'gyrase_calibration'
gyrase_k_on_min = 0.001
gyrase_k_on_max = 0.01
gyrase_k_off_min = 0.1
gyrase_k_off_max = 1.0
gyrase_k_cat_min = -20.0  # Ranges to vary k_cat
gyrase_k_cat_max = -5.0
gyrase_width_min = 0.001
gyrase_width_max = 0.05
gyrase_threshold_min = 0.001
gyrase_threshold_max = 0.05

# FIGURE
# -----------------------------------
width = 12
height = 6
lw = 3
colors = ['blue', 'red', 'green', 'black', 'grey', 'yellow', 'cyan', 'purple', 'brown', 'pink']
fig, axs = plt.subplots(2, figsize=( width, 2 * height), tight_layout=True)


# ----------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------
def Michael_Menten_equation(vmax, KM, S):
    return vmax * S / (KM + S)


def integrate_MM(vmax, KM, substrate0, product0, frames, dt):
    substrate_array = np.zeros(frames)
    product_array = np.zeros(frames)
    substrate_array[0] = substrate0
    product_array[0] = product0

    substrate = substrate0
    product = product0

    for k in range(1, frames):
        v = Michael_Menten_equation(vmax=vmax, KM=KM, S=substrate)
        product = product + v * dt
        substrate = substrate - v * dt
        substrate_array[k] = substrate
        product_array[k] = product
    return substrate_array, product_array


def rescale_product_to_sigma(product, sigma_min, sigma_max):
    #    current_min
    #    frames = len(product)
    #    sigma = np.zeros(frames)
    #    sigma = product-np.min(product)
    #    sigma = sigma/np.max(sigma) # Normalized
    #    d = sigma_f - sigma_i
    #    sigma = (sigma - sigma_i)/1#abs(sigma_f-sigma_i)

    current_min = np.min(product)
    current_max = np.max(product)
    range_ = current_max - current_min
    desired_range = sigma_max - sigma_min
    sigma = ((product - current_min) / range_) * desired_range + sigma_min
    return sigma


# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the mean superhelical
# density for each substrate concentration
def run_objective_function(list_names, list_k_cat, list_k_on, list_k_off, list_width,
                           list_threshold, list_concentration):
    # Let's run experiments for the substrate concentrations
    my_objective = 0.0
    simulation_superhelicals = []
    print(list_names, 'initial_sigma', initial_sigma)
    for s, substrate0 in enumerate(initial_substrates):
        exp_superhelical = exp_superhelicals[s]
        # Create items
        items = pt.set_items_topo_calibration(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                                              output_prefix,
                                              frames, series, continuation, dt, tm, mm, n_simulations,
                                              initial_sigma,
                                              list_names, list_k_cat, list_k_on, list_k_off, list_width, list_threshold,
                                              list_concentration, substrate0)
        # Create a multiprocessing pool
        pool = multiprocessing.Pool()
        pool_results = pool.map(pt.run_single_simulation_topo_calibration, items)

        my_supercoiling = np.zeros((frames, n_simulations))
        for i, sigma in enumerate(pool_results):
            my_supercoiling[:, i] = sigma[:-1]
        mea = np.mean(my_supercoiling, axis=1)
        my_objective += np.sum(np.square(np.mean(my_supercoiling, axis=1) - exp_superhelical))
        simulation_superhelicals.append(mea)
    my_objective = my_objective + 0.0
    return my_objective, simulation_superhelicals


def gyrase_objective_function(params):
    list_names = ['gyrase', 'topoI']

    list_k_cat, list_k_on, list_k_off, list_width, list_threshold, list_concentration = \
        set_params_lists(params)

    # Let's run experiments for the substrate concentrations
    my_objective, simulation_superhelicals = run_objective_function(list_names, list_k_cat, list_k_on, list_k_off,
                                                                    list_width, list_threshold, list_concentration)
    return my_objective


def set_params_lists(params):
    my_k_cat = params[my_vars[0]]
    my_k_on = params[my_vars[1]]
    my_k_off = params[my_vars[2]]
    my_width = params[my_vars[3]]
    my_threshold = params[my_vars[4]]

    list_k_cat = [my_k_cat, 0]
    list_k_on = [my_k_on, 0]
    list_k_off = [my_k_off, 0]
    list_width = [my_width, 0]
    list_threshold = [my_threshold, 0]
    list_concentration = [enzyme_concentration, 0]
    return list_k_cat, list_k_on, list_k_off, list_width, list_threshold, list_concentration


# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------

# Gyrase
# ==================================================================================================================
# Kinetics: RDNA + Gyrase + ATP -> RDNA-Gyrase-ATP -> SDNA + Gyrase + ADP
# Product = Less fluorescent or Supercoiled DNA
# Substrate = Concentration of Relaxed DNA or fluorescent
initial_sigma = -0.02  # Is actually the other way around, but there's an error somewhere but I'm lazy to find it
final_sigma = 0.0
initial_product = 4.0
initial_substrate = .75
initial_substrates = [0.75, 1.50, 3.6, 5.4, 7.2]
enzyme_concentration = 44.6
K_M = 2.7
k_cat = .0011
v_max = k_cat * enzyme_concentration
exp_substrates = []
exp_products = []
exp_superhelicals = []
for count, initial_substrate in enumerate(initial_substrates):
    # Substrates and products
    # ----------------------------------
    ax = axs[0]
    substrate, product = integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                      frames=frames, dt=dt)
    ax.plot(time, substrate, lw=lw, color=colors[count])

    # Sigma deduction
    # ----------------------------------
    ax = axs[1]
    superhelical = rescale_product_to_sigma(substrate, initial_sigma, final_sigma)
    ax.plot(time, superhelical, lw=lw, color=colors[count])

    exp_substrates.append(substrate)
    exp_products.append(product)
    exp_superhelicals.append(superhelical)

# Optimization
# ------------------
initial_sigma = final_sigma

space = {
    my_vars[0]: hp.uniform(my_vars[0], gyrase_k_cat_min, gyrase_k_cat_max),
    my_vars[1]: hp.uniform(my_vars[1], gyrase_k_on_min, gyrase_k_on_max),
    my_vars[2]: hp.uniform(my_vars[2], gyrase_k_off_min, gyrase_k_off_max),
    my_vars[3]: hp.uniform(my_vars[3], gyrase_width_min, gyrase_width_max),
    my_vars[4]: hp.uniform(my_vars[4], gyrase_threshold_min, gyrase_threshold_max)
}

gyrase_best = fmin(
    fn=gyrase_objective_function,  # Objective Function to optimize
    space=space,  # Hyperparameter's Search Space
    algo=tpe.suggest,  # Optimization algorithm (representative TPE)
    max_evals=tests  # Number of optimization attempts
)
print(gyrase_best)
print(gyrase_best.values())
gyrase_best_df = pd.DataFrame.from_dict(gyrase_best.values())  # TODO: fix this because it doesn't save the letters
gyrase_best_df.to_csv(gyrase_out + '.csv', index=False, sep=',')

# Plot simulation example:
# ------------------
# The idea is that the first one is the one with the concentration different to 0.
list_names = ['gyrase', 'topoI']
list_k_cat, list_k_on, list_k_off, list_width, list_threshold, list_concentration = set_params_lists(gyrase_best)
# Let's run experiments for the substrate concentrations
my_objective, simulation_superhelicals2 = run_objective_function(list_names, list_k_cat, list_k_on, list_k_off,
                                                                list_width, list_threshold, list_concentration)
# And plot
ax = axs[1]
for count, initial_substrate in enumerate(initial_substrates):
    superhelical = simulation_superhelicals2[count]
    ax.plot(time, superhelical, '--', lw=lw, color=colors[count])

for ax in [axs[0]]:
    ax.set_xlabel('time (s)')
    ax.set_ylabel('Relaxed DNA (nM)')
    ax.grid(True)

for ax in [axs[1]]:
    ax.set_xlabel('time (s)')
    ax.set_ylabel('Superhelical response')
    ax.grid(True)

axs[0].set_title('Gyrase')

plt.show()
