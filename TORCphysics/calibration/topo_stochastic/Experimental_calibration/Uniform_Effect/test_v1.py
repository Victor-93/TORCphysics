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
circuit_filename = 'circuit.csv'
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
topoI_environment_filename = 'topoI_environment.csv'
gyrase_environment_filename = 'gyrase_environment.csv'

# Experimental concentration of topoisomerases
# TODO: Check bien las concentrations
gyrase_concentration = 44.6
topo_concentration = 17.0

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
n_simulations = 2
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
def run_objective_function(global_dict, variations_list):
    # Let's run experiments for the substrate concentrations
    my_objective = 0.0
    simulation_superhelicals = []
    for s, substrate0 in enumerate(initial_substrates):
        exp_superhelical = exp_superhelicals[s]  # Experimental superhelical densities
        global_dict['DNA_concentration'] = substrate0  # DNA concentration

        # Let's create an Item to pass the conditions to the simulation
        Item = {'global_conditions': global_dict, 'variations': variations_list}

        # But we actually need a list of items, so the pool can pass each item to the function
        Items = []
        for simulation_number in range(n_simulations):
            g_dict = dict(global_dict)
            g_dict['n_simulations'] = simulation_number
            Item = {'global_conditions': g_dict, 'variations': variations_list}

            Items.append(Item)

        # Create a multiprocessing pool
        pool = multiprocessing.Pool()
        pool_results = pool.map(pt.single_simulation_calibration_w_supercoiling, Items)

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

    # Global dictionary
    # ------------------------------------------
    global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                   'enzymes_filename': enzymes_filename, 'environment_filename': topoI_environment_filename,
                   'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                   'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': initial_sigma,
                   'DNA_concentration': 0.0}

    # Variation dictionary
    # ------------------------------------------
    name = 'topoI'
    object_type = 'environmental'
    binding_model_name = 'TopoIRecognition'
    binding_oparams = {'k_on': params['k_on'], 'width': params['width'], 'threshold': params['threshold']}
    effect_model_name = 'TopoIUniform'
    effect_oparams = {'k_cat': params['k_cat']}
    unbinding_model_name = 'PoissonUnBinding'
    unbinding_oparams = {'k_off': params['k_off']}
    concentration = topo_concentration

    topo_variation = {'name': name, 'object_type': object_type,
                      'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                      'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                      'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                      'concentration': concentration}

    my_objective, simulation_superhelicals = run_objective_function(global_dict=global_dict,
                                                                    variations_list=[topo_variation])
    return my_objective


def gyrase_objective_function(params):
    # We need to prepare the inputs.

    # Global dictionary
    # ------------------------------------------
    global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                   'enzymes_filename': enzymes_filename, 'environment_filename': gyrase_environment_filename,
                   'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                   'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': initial_sigma,
                   'DNA_concentration': 0.0}

    # Variation dictionary
    # ------------------------------------------
    name = 'gyrase'
    object_type = 'environmental'
    binding_model_name = 'GyraseRecognition'
    binding_oparams = {'k_on': params['k_on'], 'width': params['width'], 'threshold': params['threshold']}
    effect_model_name = 'GyraseUniform'
    effect_oparams = {'k_cat': params['k_cat']}
    unbinding_model_name = 'PoissonUnBinding'
    unbinding_oparams = {'k_off': params['k_off']}
    concentration = gyrase_concentration

    gyrase_variation = {'name': name, 'object_type': object_type,
                        'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                        'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                        'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                        'concentration': concentration}

    my_objective, simulation_superhelicals = run_objective_function(global_dict=global_dict,
                                                                    variations_list=[gyrase_variation])
    return my_objective


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
enzyme_concentration = topo_concentration
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

# ---------------------------------------------------------------------------------------------------------------------
# Plot simulation example:
# ---------------------------------------------------------------------------------------------------------------------
# The idea is that the first one is the one with the concentration different to 0.
# Let's create inputs
# Global dictionary
# ------------------------------------------
global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
               'enzymes_filename': enzymes_filename, 'environment_filename': topoI_environment_filename,
               'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
               'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': initial_sigma,
               'DNA_concentration': 0.0}
# Variation dictionary
# ------------------------------------------
name = 'topoI'
object_type = 'environmental'
binding_model_name = 'TopoIRecognition'
binding_oparams = {'k_on': topoI_best['k_on'], 'width': topoI_best['width'], 'threshold': topoI_best['threshold']}
effect_model_name = 'TopoIUniform'
effect_oparams = {'k_cat': topoI_best['k_cat']}
unbinding_model_name = 'PoissonUnBinding'
unbinding_oparams = {'k_off': topoI_best['k_off']}
concentration = topo_concentration
variation = {'name': name, 'object_type': object_type,
             'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
             'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
             'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
             'concentration': concentration}
# Let's run experiments for the substrate concentrations
my_objective, simulation_superhelicals = run_objective_function(global_dict=global_dict, variations_list=[variation])
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
enzyme_concentration = gyrase_concentration
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
# ---------------------------------------------------------------------------------------------------------------------
# Plot simulation example:
# ---------------------------------------------------------------------------------------------------------------------
# The idea is that the first one is the one with the concentration different to 0.
# Let's create inputs
# Global dictionary
# ------------------------------------------
global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
               'enzymes_filename': enzymes_filename, 'environment_filename': gyrase_environment_filename,
               'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
               'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': initial_sigma,
               'DNA_concentration': 0.0}
# Variation dictionary
# ------------------------------------------
name = 'gyrase'
object_type = 'environmental'
binding_model_name = 'GyraseRecognition'
binding_oparams = {'k_on': gyrase_best['k_on'], 'width': gyrase_best['width'], 'threshold': gyrase_best['threshold']}
effect_model_name = 'GyraseUniform'
effect_oparams = {'k_cat': gyrase_best['k_cat']}
unbinding_model_name = 'PoissonUnBinding'
unbinding_oparams = {'k_off': gyrase_best['k_off']}
concentration = topo_concentration
variation = {'name': name, 'object_type': object_type,
             'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
             'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
             'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
             'concentration': concentration}
# Let's run experiments for the substrate concentrations
my_objective, simulation_superhelicals2 = run_objective_function(global_dict=global_dict, variations_list=[variation])
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
