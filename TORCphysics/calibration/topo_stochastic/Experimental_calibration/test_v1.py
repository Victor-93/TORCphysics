import matplotlib.pyplot as plt
import numpy as np
from hyperopt import tpe, hp, fmin
import multiprocessing
from TORCphysics import parallelization_tools as pt
import pandas as pd
import topo_calibration_tools as tct

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
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
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
environment_filename = 'environment_test.csv'
# Antiguo environment_test:
# type,name,site_type,concentration,k_on,k_off,k_cat,size,oparams
# topo,topoI,DNA,0.25,0.005,0.5,7.5,120,none
# topo,gyrase,DNA,0.25,0.005,0.5,-12.5,120,none


tm = 'stochastic'
output_prefix = 'test0'
series = True
continuation = False
mm = 'uniform'

# For parallelization and calibration
n_simulations = 1
tests = 2  # 00  # number of tests for parametrization

my_vars = ['k_cat', 'k_on', 'k_off', 'width', 'threshold']

# RANGES FOR RANDOM SEARCH
# -----------------------------------
# Gyrase ranges
gyrase_out = 'gyrase_calibration_test'
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

# TopoI ranges
topoI_out = 'topoI_calibration_test'
topoI_k_on_min = 0.001
topoI_k_on_max = 0.01
topoI_k_off_min = 0.1
topoI_k_off_max = 1.0
topoI_k_cat_min = 5.0  # Ranges to vary k_cat
topoI_k_cat_max = 20.0
topoI_width_min = 0.001
topoI_width_max = 0.05
topoI_threshold_min = -0.05
topoI_threshold_max = -0.001

# Meyer topos curve parameters:
topo_w = 0.012  # width
topo_t = -0.04  # thresholds
topo_k = 0.001  # k_cat
gyra_w = 0.025  # width
gyra_t = 0.01  # threshold
gyra_k = 0.001  # k_cat

# FIGURE
# -----------------------------------
width = 8
height = 4
lw = 3
colors = ['blue', 'red', 'green', 'black', 'grey', 'yellow', 'cyan', 'purple', 'brown', 'pink']
fig, axs = plt.subplots(2, 2, figsize=(2 * width, 2 * height), tight_layout=True)


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
        global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                       'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                       'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                       'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': initial_sigma,
                       'DNA_concentration': substrate0}

        # Create items
        items = pt.set_items_topo_calibration(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                                              output_prefix,
                                              frames, series, continuation, dt, n_simulations,
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


def topoI_objective_function(params):
    #   * So what we have is a dictionary with many entries.
    #   'initial_conditions': dict with global parameters, including files (environment, enzymes, sites, etc...)
    #   * We also have 'variations': list with dictionaries. Every entry is regarding a molecule/model to calibrate.
    #   * In general, the idea is that an external costumable function builds this type of entry. So we have.
    #   Item = {'initial_conditions': dict{superhelical:-0.04, 'circuit_filename': csv, 'dt':0.1, etc..}.
    #           'variations': list[ dict{'name': topoI, 'type': environment, 'binding_model': TopoI,
    #           'effect_model': None, etc...}

    # We need to prepare the inputs.
    name = 'topoI'
    mol_type = 'environmental'
    binding_model_name = 'TopoIRecognition'
    binding_oparams = {'k_on': params['k_on'], 'width': params['width'], 'threshold': params['threshold']}
    effect_model_name = 'TopoIUniform'
    effect_oparams = {'k_cat': params['k_cat']}
    unbinding_model_name = 'PoissonUnBinding'
    unbinding_oparams = {'k_off': params['k_off']}
    concentration = 0.2

    topo_variation = {'name': name, 'mol_type': mol_type,
                      'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                      'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                      'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                      'concentration': concentration}
    # TODO: My idea here was to brindg the objective function here? Or maybe just find a way to fill th econcentrations.
    #  I think it doesn't make a difference if I bring the objective function here, it might be easier that way.

    # TODO: IDEA. Notice that, I could also try to define models here, in this script... If I want to define new
    #  models, I could just do it here... Load it once, and then just have the same model for every simulation test.
    #  Or maybe it doesn't make that much difference. You know what, I'll leave it just as it is now.
    # Let's dp the variations list of dicts
    variations_list = []
    # Let's start with the globals

    #    circuit_filename, sites_filename, enzymes_filename, environment_filename,
    #    output_prefix,
    ##    frames, series, continuation, dt, n_simulations,
    #    initial_sigma,
    #    list_names, list_k_cat, list_k_on, list_k_off, list_width, list_threshold,
    #    list_concentration, substrate0)

    list_names = ['topoI', 'gyrase']

    list_k_cat, list_k_on, list_k_off, list_width, list_threshold, list_concentration = \
        set_params_lists(params)

    # Let's run experiments for the substrate concentrations
    my_objective, simulation_superhelicals = run_objective_function(list_names, list_k_cat, list_k_on, list_k_off,
                                                                    list_width, list_threshold, list_concentration)
    return my_objective


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

# TOPOISOMERASE I
# ==================================================================================================================
# Kinetics: SDNA + TopoI -> SDNA-TopoI -> RDNA + TopoI
# Product = Fluorescent or Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs
initial_sigma = -.047
final_sigma = 0.0
initial_product = 0.0
initial_substrate = .7
initial_substrates = [0.35, 0.7, 1.1, 1.8, 2.1]
enzyme_concentration = 17
K_M = 1.5
k_cat = .0023  # 0.003
v_max = k_cat * enzyme_concentration
exp_substrates = []
exp_products = []
exp_superhelicals = []
for count, initial_substrate in enumerate(initial_substrates):
    # Substrates and products
    # ----------------------------------
    ax = axs[0, 0]
    substrate, product = tct.integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                          frames=frames, dt=dt)
    ax.plot(time, product, lw=lw, color=colors[count])

    # Sigma deduction
    # ----------------------------------
    ax = axs[1, 0]
    superhelical = tct.rescale_product_to_sigma(product, initial_sigma, final_sigma)
    ax.plot(time, superhelical, lw=lw, color=colors[count])

    exp_substrates.append(substrate)
    exp_products.append(product)
    exp_superhelicals.append(superhelical)

# TODO: Here is what you need to fix:
#  1.- When you do the space thing and call your function, you already have params; so you can give this dictionary
#   directly to your binding, effect or unbinding model.
#  2.- So we have to pass those params to the model.oparams.
#   2.1.- We also need to give some initial conditions to the simulations.
#   * So what we have is a dictionary with many entries.
#   'initial_conditions': dict with global parameters, including files (environment, enzymes, sites, etc...)
#   * We also have 'variations': list with dictionaries. Every entry is regarding a molecule/model to calibrate.
#   * In general, the idea is that an external costumable function builds this type of entry. So we have.
#   Item = {'initial_conditions': dict{superhelical:-0.04, 'circuit_filename': csv, 'dt':0.1, etc..}.
#           'variations': list[ dict{'name': topoI, 'type': environment, 'binding_model': TopoI,
#           'effect_model': None, etc...}
#           }
#  3.- We pass it to a function in parallelization that sets up all of this.
# Optimization
# ------------------
space = {
    my_vars[0]: hp.uniform(my_vars[0], topoI_k_cat_min, topoI_k_cat_max),
    my_vars[1]: hp.uniform(my_vars[1], topoI_k_on_min, topoI_k_on_max),
    my_vars[2]: hp.uniform(my_vars[2], topoI_k_off_min, topoI_k_off_max),
    my_vars[3]: hp.uniform(my_vars[3], topoI_width_min, topoI_width_max),
    my_vars[4]: hp.uniform(my_vars[4], topoI_threshold_min, topoI_threshold_max)
}

topoI_best = fmin(
    fn=topoI_objective_function,  # Objective Function to optimize
    space=space,  # Hyperparameter's Search Space
    algo=tpe.suggest,  # Optimization algorithm (representative TPE)
    max_evals=tests  # Number of optimization attempts
)
print(topoI_best)
print(topoI_best.values())
topoI_best_df = pd.DataFrame.from_dict(topoI_best.values())  # TODO: fix this because it doesn't save the letters
topoI_best_df.to_csv(topoI_out + '.csv', index=False, sep=',')
## Read dictionary file (should be .pkl)
# with open(topoI_out, 'wb') as fp:
#    pickle.dump(topoI_best, fp)
#    print('dictionary saved successfully to file')

# Plot simulation example:
# ------------------
# The idea is that the first one is the one with the concentration different to 0.
list_names = ['topoI', 'gyrase']
list_k_cat, list_k_on, list_k_off, list_width, list_threshold, list_concentration = set_params_lists(topoI_best)
# Let's run experiments for the substrate concentrations
my_objective, simulation_superhelicals = run_objective_function(list_names, list_k_cat, list_k_on, list_k_off,
                                                                list_width, list_threshold, list_concentration)
# And plot
ax = axs[1, 0]
for count, initial_substrate in enumerate(initial_substrates):
    superhelical = simulation_superhelicals[count]
    ax.plot(time, superhelical, '--', lw=lw, color=colors[count])

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
    ax = axs[0, 1]
    substrate, product = tct.integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                          frames=frames, dt=dt)
    ax.plot(time, substrate, lw=lw, color=colors[count])

    # Sigma deduction
    # ----------------------------------
    ax = axs[1, 1]
    superhelical = tct.rescale_product_to_sigma(substrate, initial_sigma, final_sigma)
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
ax = axs[1, 1]
for count, initial_substrate in enumerate(initial_substrates):
    superhelical = simulation_superhelicals2[count]
    ax.plot(time, superhelical, '--', lw=lw, color=colors[count])

for ax in [axs[0, 0], axs[0, 1]]:
    ax.set_xlabel('time (s)')
    ax.set_ylabel('Relaxed DNA (nM)')
    ax.grid(True)

for ax in [axs[1, 0], axs[1, 1]]:
    ax.set_xlabel('time (s)')
    ax.set_ylabel('Superhelical response')
    ax.grid(True)

axs[0, 0].set_title('Topoisomerase I')
axs[0, 1].set_title('Gyrase')

plt.savefig('calibration_test.png')
plt.show()
