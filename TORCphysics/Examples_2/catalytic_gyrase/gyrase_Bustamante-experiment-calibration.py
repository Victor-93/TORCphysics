import sys
#sys.path.append("/mnt/parscratch/users/username")  # For stanage #TODO: Modify this to your directory, where you uploaded TORCphysics
import numpy as np
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em
from TORCphysics import unbinding_model as unbm
import parallelization_module as pm
from custom_gyrase import GyraseCyclesEffect, GyraseCyclesUnbinding, GyraseCyclesForce
import objective_single_molecule_experiment as osme
from TORCphysics import Site
import pickle
import copy

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This script performs a parameter sweep to find the best parameterisation set that reproduces the observations from
# Bustamante et al: Mechanochemical analysis of DNA gyrase using rotor bead tracking.

# TIP: Check the initial conditions, and make sure that they make sense and that you request the
# necessary HPC resources.
# Also check the parameter space so the parameters explored make sense.

# The output of this file is the hyperopt object (*trials.pkl), which contains all the information
# for parameter tests and their results.
# Check the analysis script with the debug option so you can see for yourself the shape of the data.
# Also, check the output of the objective function so you get a sense of the information stored,
# which essentially is the loss (error between experiment and simulation, and the rotations induced
# per binding event for each force tested.

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
file_out = 'gyrase_calibration' #+ '_' + str(dt)
tests = 10000  # 10000   # number of tests for parametrization
#parallel = False
parallel = True
n_simulations = 50 #2  # 100 # For stanage

# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 0.1
#dt = 1.0

# Time parameters
initial_time = 0

# For the second experiment
final_time2 = 10800
#final_time2 = 10800 // 5
time2 = np.arange(initial_time, final_time2 + dt, dt)
frames2 = len(time2)

# Concentrations in nM
gyrase_concentration = 10.0
#gyrase_concentration = 44.6

# Superhelical values (sigma) for each case
sigma0_RX = 0.0  # Superhelical density of relaxed DNA
# At this value the torque is too strong.

# Models to calibrate to calibrate
# -----------------------------------
gyrase_binding_model = bm.PoissonBinding()
gyrase_effect_model = GyraseCyclesForce()
gyrase_unbinding_model = GyraseCyclesUnbinding()

gyrase_k_cat = 20.0 # 21 # How many turns are added per catalitic event

single_molecule_size = 2200
Forces = [0.35, 0.6, 0.8, 1.0, 1.3] # pN
#Forces = [0.35, 1.3] # pN


# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the average
# number of rotations induced per binding event. This can be used to quantify the number of cycles performed

def objective_function(params):
    # Sort parameters and load them to the models
    # ------------------------------------------------------------------
    objective = 0
    total_rotations = []
    for fi, my_force in enumerate(Forces):  # fi is the index
        gyrase_params = pm.extract_params(params, 'gyrase')
        gyrase_params['force'] = my_force
        gyrase_params['k_cat'] = gyrase_k_cat

        gyrase_bm = gyrase_binding_model.__class__(**gyrase_params)
        gyrase_em = gyrase_effect_model.__class__(**gyrase_params)
        gyrase_unbm = gyrase_unbinding_model.__class__(**gyrase_params)

        # -----------------------------------------------------------------------------------------------
        # Experiment
        # -----------------------------------------------------------------------------------------------
        # This experiment is based on the "Mechanochemical analysis of DNA gyrase using rotor bead tracking" paper.
        # We model gyrase binding only to one binding site.
        # Here, we test two systems:
        # 1.- Gyrase acting on relaxed DNA
        # 2.- Topoisomerase I acting on negatively supercoiled DNA (we don't have data for this one).

        # In this experiment at the moment we don't collect errors, but we might in the future since in principle,
        # we could fit Gaussians

        # Define the sites for both enzymes
        # ------------------------------------------
        gyrase_site = Site(site_type='gyrase', name='gyrase_s', start=1000, end=0,
                           binding_model=copy.deepcopy(gyrase_bm))

        # Global dictionaries - One for each system - just one system for now to keep it simple
        # ------------------------------------------
        global_dict_gyrase2 = {'frames': frames2, 'dt': dt, 'superhelicity': sigma0_RX, 'system_name': 'gyrase_RX',
                               'size': single_molecule_size, 'structure': 'linear',
                               'gyrase_concentration': gyrase_concentration,
                               'topoI': False, 'gyrase': True,
                               'gyrase_bm': gyrase_bm, 'gyrase_em': gyrase_em, 'gyrase_unbm': gyrase_unbm,
                               'sites': [gyrase_site]
                               }
        # Run experiments
        # ----------------------------------------------------
        global_dict_list2 = [global_dict_gyrase2]

        exp2_output = pm.run_experiment_2(global_dict_list=global_dict_list2, n_simulations=n_simulations,
                                          parallelization=parallel)

        # Here calculate the objective of the current force
        sim_objective, sim_rotations = osme.Bustamante_objective(exp2_output, my_force)

        # Add it to the overall objective
        total_rotations.append(sim_rotations)
        objective += sim_objective

    # Arrange outputs
    output_dict = {'loss': objective, 'status': STATUS_OK,  # hyperopt results
                   # Specific results for processing
                   'rotations_distribution': total_rotations}
    # loss: is the sum of squared errors between experimental curve and simulation at all forces
    # rotations_distribution: It is a list that contains the rotation distributions. Each entry represents
    #                         the rotations at a particular force, so remember the forces tested.
    return output_dict

#----------------------------------------------------------------------------------------------------------------------
# Random Search algorithm!
#----------------------------------------------------------------------------------------------------------------------
trials = Trials() # We need this for hyperopt!

# Here, you define the parameter space for the random search. These are basically the limits and parameter sampling.
# For simplicity, we use a uniform distribution for sampling each parameter, and hopefully we explore the "complete"
# parameter space (which is multidimensional).
space = {

    # Gyrase params
    # ----------------------------------
    # Binding
    'k_on_gyrase': hp.uniform('k_on_gyrase', 0.001, 0.1),
    #'k_on_gyrase': hp.uniform('k_on_gyrase', 0.000024, 0.00024),
    # 'width_gyrase': hp.uniform('width_gyrase', 0.00125, 0.0125),
    # 'threshold_gyrase': hp.uniform('threshold_gyrase', -0.06, 0.04),

    # Effect
    # 'k_cat_gyrase': hp.uniform('k_cat_gyrase', 0.0, 20.0),
    # We don't have a sense of these rates...
    'k_wrap_gyrase': hp.uniform('k_wrap_gyrase', 0.01, 1.0),
    'k_unwrap_gyrase': hp.uniform('k_unwrap_gyrase', 0.01, 1.0),
    'k_go_gyrase': hp.uniform('k_go_gyrase', 0.01, 1.0),
    'x_wrap_gyrase': hp.uniform('x_wrap_gyrase', 0.01, 30.0),

    # Unbinding
    'k_off_gyrase': hp.uniform('k_off_gyrase', 0.01, 1.0)
}

# Define the file where you want to save the output
output_file_path = file_out + '.info'

# Open the file in write mode
with open(output_file_path, 'w') as f:
    # Redirect the standard output to the file
    sys.stdout = f  # Comment this one out to print stuff to your screen.

    # Your code that prints to the screen
    print("Hello, this is the info file for the calibration of Topo I and Gyrase Models.")
    print("Launching calibration for dt="+str(dt))
    print('Ran ' + str(n_simulations) + ' simulations per test. ')
    print('Number of tests = ' + str(tests))

    best = fmin(
        fn=objective_function,  # Objective Function to optimize
        space=space,  # Hyperparameter's Search Space
        algo=tpe.suggest,  # Optimization algorithm (representative TPE)
        max_evals=tests,  # Number of optimization attempts
        trials=trials
    )

    print(" ")
    print("Optimal parameters found from random search: ")
    print(best)

# -------------------------------------------------------------------------------------------
# Output the useful information
# -------------------------------------------------------------------------------------------

# Save the trials object! - this object contains all the information
# --------------------------------------------------------------------------
with open(file_out+'-trials.pkl', 'wb') as file:
    pickle.dump(trials, file)

# Save best parameterisations as csv
# --------------------------------------------------------------------------
# If you want to save all best params in a single csv ...
#best_df = pd.DataFrame.from_dict([best])
#best_df.to_csv(file_out + '.csv', index=False, sep=',')

# Let's save the best parameterisation for each enzyme on a separate file each
for name in ['gyrase']:
    my_params = pm.extract_params(best, name)
    my_df = pd.DataFrame.from_dict([my_params])
    my_df.to_csv(file_out+'best_'+name+'.csv', index=False, sep=',')

# Let's extract the losses related to each parameterisation, and save it as csv to perform quick statistics.
# --------------------------------------------------------------------------
keys = list(best.keys())
for n in range(tests):
    tdi = trials.trials[n]  # dictionary with results for test n
    va = trials.trials[n]['misc']['vals']  #values
    my_dict = {}
    my_dict['test'] = n
    my_dict['loss'] = trials.trials[n]['result']['loss']  # loss
    for key in keys:
        my_dict[key] = va[key][0]

    # Add a new row using append method
    if n == 0:
        params_df = pd.DataFrame.from_dict([my_dict])
    else:
        new_row = pd.DataFrame.from_dict([my_dict])
        params_df = pd.concat([params_df, new_row], ignore_index=True)

params_df.to_csv(file_out+'_values.csv', index=False, sep=',')

