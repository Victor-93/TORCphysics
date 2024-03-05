import numpy as np
from hyperopt import tpe, hp, fmin
import pandas as pd
import topo_calibration_tools as tct
import sys

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
environment_filename = 'gyrase_environment.csv'

# Experimental concentration of topoisomerases
# gyrase_concentration = 44.6
mol_concentration = 44.6

# Relaxed is the superhelical density when the DNA is relaxed
relaxed_DNA = 0.01
final_DNA = -0.1  # This is after the enzyme effect

tm = 'stochastic'
output_prefix = 'test0'
series = True
continuation = False
mm = 'uniform'

# For parallelization and calibration
n_simulations = 60 #48 #120
tests = 100  # number of tests for parametrization

# Molecule/model to calibrate
# -----------------------------------
mol_name = 'gyrase'
mol_type = 'environmental'
mol_binding_model_name = 'PoissonBinding'
mol_effect_model_name = 'GyraseLinear'
mol_unbinding_model_name = 'PoissonUnBinding'

# RANGES FOR RANDOM SEARCH
# -----------------------------------
# Gyrase ranges
file_out = mol_name + '_calibration'
k_on_min = 0.001
k_on_max = 0.1
k_off_min = 0.01
k_off_max = 1.0
k_cat_min = 5.0  # Ranges to vary k_cat
k_cat_max = 20.0
sigma0_min = -0.3
sigma0_max = 0.0

# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the mean superhelical
# density for each substrate concentration


def objective_function(params):
    # We need to prepare the inputs.

    # Global dictionary
    # ------------------------------------------
    global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                   'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                   'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                   'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': initial_sigma,
                   'DNA_concentration': 0.0}

    # Variation dictionary
    # ------------------------------------------
    name = mol_name
    object_type = mol_type
    binding_model_name = mol_binding_model_name
    binding_oparams = {'k_on': params['k_on']}
    effect_model_name = mol_effect_model_name
    effect_oparams = {'k_cat': params['k_cat'], 'sigma0': params['sigma0']}
    unbinding_model_name = mol_unbinding_model_name
    unbinding_oparams = {'k_off': params['k_off']}
    concentration = mol_concentration  # / mol_concentration  # Because this is the reference.

    topo_variation = {'name': name, 'object_type': object_type,
                      'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                      'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                      'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                      'concentration': concentration}

    my_objective, simulation_superhelicals = tct.run_objective_function(global_dict=global_dict,
                                                                        variations_list=[topo_variation],
                                                                        initial_substrates=initial_substrates,
                                                                        exp_superhelicals=exp_superhelicals,
                                                                        n_simulations=n_simulations)
    return my_objective


# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------

# Experimental evaluation
# ==================================================================================================================
# Kinetics: SDNA + TopoI -> SDNA-TopoI -> RDNA + TopoI
# Product = Fluorescent or Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs
# Kinetics: SDNA + TopoI -> SDNA-TopoI -> RDNA + TopoI
# Product = Fluorescent or Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs
initial_sigma = final_DNA  # Is actually the other way around, but there's an error somewhere but I'm lazy to find it
final_sigma = relaxed_DNA  # 0.0
initial_product = 4.0
initial_substrate = .75
initial_substrates = [0.75]
enzyme_concentration = mol_concentration
K_M = 2.7
k_cat = .0011
v_max = k_cat * enzyme_concentration
exp_substrates = []
exp_products = []
exp_superhelicals = []
for count, initial_substrate in enumerate(initial_substrates):
    # Substrates and products
    # ----------------------------------
    substrate, product = tct.integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                          frames=frames, dt=dt)

    # Sigma deduction
    # ----------------------------------
    superhelical = tct.rescale_product_to_sigma(substrate, initial_sigma, final_sigma)

    # Collect results
    exp_substrates.append(substrate)
    exp_products.append(product)
    exp_superhelicals.append(superhelical)

# Optimization
# ==================================================================================================================
initial_sigma = final_sigma

space = {
    'k_cat': hp.uniform('k_cat', k_cat_min, k_cat_max),
    'k_on': hp.uniform('k_on', k_on_min, k_on_max),
    'k_off': hp.uniform('k_off', k_off_min, k_off_max),
    'sigma0': hp.uniform('sigma0', sigma0_min, sigma0_max)
}

# Save the current standard output
original_stdout = sys.stdout
# Define the file where you want to save the output
output_file_path = file_out + '.info'

# Open the file in write mode
with open(output_file_path, 'w') as f:
    # Redirect the standard output to the file
    sys.stdout = f

    # Your code that prints to the screen
    print("Hello, this is the info file for the calibration of " + mol_name + ' model.')
    print("Binding Model = " + mol_binding_model_name)
    print("Effect Model = " + mol_effect_model_name)
    print("Unbinding Model = " + mol_unbinding_model_name)
    print('Ran ' + str(n_simulations) + ' simulations per test. ')
    print('Number of tests = ' + str(tests))

    best = fmin(
        fn=objective_function,  # Objective Function to optimize
        space=space,  # Hyperparameter's Search Space
        algo=tpe.suggest,  # Optimization algorithm (representative TPE)
        max_evals=tests  # Number of optimization attempts
    )

    print(" ")
    print("Optimal parameters found from random search: ")
    print(best)

best_df = pd.DataFrame.from_dict([best])  # TODO: fix this because it doesn't save the letters
best_df.to_csv(file_out + '.csv', index=False, sep=',')
