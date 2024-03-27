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

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 0.25
initial_time = 0
final_time = 600
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'calibration'

# For the simulation
circuit_filename = 'circuit.csv'
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
environment_filename = 'topoI_environment.csv'
gyrase_environment_filename = 'gyrase_environment.csv'

# Concentrations in nM
DNA_concentration = 0.75
gyrase_concentration = 44.6
topoI_concentration = 17.0

# MM kinetics
K_M_topoI = 1.5
k_cat_topoI = .0023
v_max_topoI = k_cat_topoI * topoI_concentration
K_M_gyrase = 2.7
k_cat_gyrase = .0011
v_max_gyrase = k_cat_gyrase * gyrase_concentration

# Superhelical values (sigma) for each case
sigma_0_topo = -0.075  # Approximately -20 supercoils according the paper
sigma_0_gyrase = 0.0  # We suppose this one.
sigma_f_gyrase = -0.1  # We also assume this one, which is the maximum at which gyrase acts.
# At this value the torque is too strong.

output_prefix = 'test0'
series = True
continuation = False

# For parallelization and calibration
n_simulations = 6#60 #48 #120
tests = 5#100  # number of tests for parametrization

# Molecule/model to calibrate
# -----------------------------------
mol_name = 'topoI'
mol_type = 'environmental'
mol_binding_model_name = 'TopoIRecognition'
mol_effect_model_name = 'TopoILinear'
mol_unbinding_model_name = 'PoissonUnBinding'
mol_concentration = topoI_concentration

# RANGES FOR RANDOM SEARCH
# -----------------------------------
# TopoI ranges
file_out = mol_name + '_calibration'
k_on_min = 0.0001
k_on_max = 0.01
k_off_min = 0.01
k_off_max = 1.0
k_cat_min = 5.0  # Ranges to vary k_cat
k_cat_max = 20.0
width_min = 0.001
width_max = 0.05
threshold_min = -0.05
threshold_max = -0.001

# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the mean superhelical
# density for each substrate concentration

def objective_function_topoI(params):
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
    binding_oparams = {'k_on': params['k_on'], 'width': params['width'], 'threshold': params['threshold']}
    effect_model_name = mol_effect_model_name
    effect_oparams = {'k_cat': params['k_cat']}
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
                                                                        exp_superhelicals=exp_sigma,
                                                                        n_simulations=n_simulations)
    return my_objective


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# ==================================================================================================================
# ==================================================================================================================
# TOPO I !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ==================================================================================================================
# ==================================================================================================================

# Build experimental curve for TOPO I
# -----------------------------------------------------
# Kinetics: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI
# Product = Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs; which initially is the same as the DNA concentration

# Integrate MM kinetics
# ------------------------------------------
# Initially, there's no relaxed DNA, and all the supercoiled DNA concentration corresponds to the plasmid conc.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_topoI(vmax=v_max_topoI, KM=K_M_topoI,
                                                      Supercoiled_0=DNA_concentration, Relaxed_0=0.0,
                                                      frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
exp_sigma = tct.topoI_to_sigma(Relaxed=relaxed_DNA, DNA_concentration=DNA_concentration, sigma0=sigma_0_topo)


# TODO: This is wrong, the objective function needs to consider the three experiments, topo I alone, Gyrase alone,
#       and both enzymes acting together.
# TODO: Calculate exp_sigma for the three cases
# TODO: Create optimization space that varies both topo I and Gyrase params
# TODO:  Create objective function that sums ob_topo + ob_gyrase + ob_both
# Optimization
# -----------------------------------------------------
space = {
    'k_cat': hp.uniform('k_cat', k_cat_min, k_cat_max),
    'k_on': hp.uniform('k_on', k_on_min, k_on_max),
    'k_off': hp.uniform('k_off', k_off_min, k_off_max),
    'width': hp.uniform('width', width_min, width_max),
    'threshold': hp.uniform('threshold', threshold_min, threshold_max)
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
