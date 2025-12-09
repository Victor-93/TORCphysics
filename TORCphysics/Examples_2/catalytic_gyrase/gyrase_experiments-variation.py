import sys
import numpy as np
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials
import pandas as pd
import pickle
import copy
import parallelization_module as pm

#sys.path.append("/mnt/parscratch/users/username")  # For stanage
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em
from TORCphysics import unbinding_model as unbm
from custom_gyrase import GyraseCyclesEffect, GyraseCyclesUnbinding
from TORCphysics import Site

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# The idea is to run the topoisomerases experiment described in the TORCphysics paper but using the new gyrase model.
# The experiment is based on the paper:
#   Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching
# The idea now, is that we will use the kinetics models from that experiment to define an error function, that will
# help us to evaluate the gyrase model performance.
# Ideally, the new model should be able to reproduce the results from that Kinetic experiment, as well as the
# single-molecule data from Bustamante et al: Mechanochemical analysis of DNA gyrase using rotor bead tracking

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0
tests = 2 #10000   # number of tests for parametrization
parallel = False
#parallel = True  # This option doesn't work on windows

# For parallelization and calibration
n_simulations = 2 #100 # For stanage

# Time parameters
initial_time = 0

# For the first experiment
final_time = 50#0  #500 # 600
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# For the second experiment
#final_time2 = 10800
final_time2 = 10800 // 10
time2 = np.arange(initial_time, final_time2 + dt, dt)
frames2 = len(time2)

file_out = 'topoisomerases_calibration' #+ '_' + str(dt)

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
sigma0_RX = 0.0  # Superhelical density of relaxed DNA
sigma0_NSC = -0.11  # Superhelical density of negatively supercoiled DNA
sigma_0_topo = -0.11  #-0.076  # Approximately -20 supercoils according the paper
sigma_0_gyrase = 0.0  # We suppose this one.
sigma_f_gyrase = -0.11  # We also assume this one, which is the maximum at which gyrase acts.
# At this value the torque is too strong.

output_prefix = 'test0'
series = True
continuation = False

# Models to calibrate to calibrate
# -----------------------------------
topoI_binding_model = bm.TopoIRecognition()  # Binding Model
topoI_effect_model = em.TopoILinear() # Effect Model
topoI_unbinding_model = unbm.PoissonUnBinding() # Unbinding Model

gyrase_binding_model = bm.GyraseRecognition()
gyrase_effect_model = GyraseCyclesEffect()
gyrase_unbinding_model = GyraseCyclesUnbinding()

# Sizes
topoI_size = 20
gyrase_size = 30

# For experiment 1
plasmid_size = 2757
topoI_ngrid = int(plasmid_size/topoI_size)
gyrase_ngrid = int(plasmid_size/gyrase_size)

# For experiment 2
single_molecule_size = 2200

# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the mean superhelical
# density for each substrate concentration

def objective_function(params):

    # -----------------------------------------------------------------------------------------------
    # First, let's test the k_on/k_off ratio.
    # This idea is to filter parameters for which the ratio is higher than one.
    # Basically, we don't want more than 1 topoisomerase binding per second (preferably, as this is random)
    # ------------------------------------------------------------------
    test_topoI_kon = (params['k_on_topoI'] * topoI_ngrid * topoI_concentration)
    test_gyrase_kon =  (params['k_on_gyrase'] * gyrase_ngrid * gyrase_concentration)

    topoI_ratio = test_topoI_kon / params['k_off_topoI']
    gyrase_ratio = test_gyrase_kon / params['k_off_gyrase']

    # Test ratio - if it's bigger than one, then we don't want it
    #if topoI_ratio > 1.0 or gyrase_ratio > 1.0:
    #    return 100 # Big objective

    # Sort parameters and load them to the models
    # ------------------------------------------------------------------
    topoI_params = pm.extract_params(params, 'topoI')
    gyrase_params = pm.extract_params(params, 'gyrase')

    topoI_bm = topoI_binding_model.__class__(**topoI_params)
    topoI_em = topoI_effect_model.__class__(**topoI_params)
    topoI_unbm = topoI_unbinding_model.__class__(**topoI_params)

    gyrase_bm = gyrase_binding_model.__class__(**gyrase_params)
    gyrase_em = gyrase_effect_model.__class__(**gyrase_params)
    gyrase_unbm = gyrase_unbinding_model.__class__(**gyrase_params)

    # -----------------------------------------------------------------------------------------------
    # Experiment 1!
    # -----------------------------------------------------------------------------------------------
    # Here, we want to calculate an error to evaluate each parameterisation set.
    # We test four different systems:
    # 1.- Topoisomerase I acting on supercoiled DNA
    # 2.- Gyrase acting on Relaxed DNA
    # 3.- Both enzymes acting on supercoiled
    # 4.- Both enzymes acting on relaxed DNA

    # Global dictionaries - One for each system
    # ------------------------------------------
    global_dict_topoI = {'frames': frames, 'dt': dt, 'superhelicity': sigma_0_topo, 'system_name': 'topoI',
                         'size': plasmid_size, 'structure': 'circular',
                         'topoI_concentration': topoI_concentration,
                         'topoI': True, 'gyrase': False, # If we want to include topoI and gyrase
                         'topoI_bm': topoI_bm, 'topoI_em': topoI_em, 'topoI_unbm': topoI_unbm
                         }

    global_dict_gyrase = {'frames': frames, 'dt': dt, 'superhelicity': sigma_0_gyrase, 'system_name': 'gyrase',
                         'size': plasmid_size, 'structure': 'circular',
                         'gyrase_concentration': gyrase_concentration,
                         'topoI': False, 'gyrase': True,
                         'gyrase_bm': gyrase_bm, 'gyrase_em': gyrase_em, 'gyrase_unbm': gyrase_unbm
                         }

    global_dict_both_sc = {'frames': frames, 'dt': dt, 'superhelicity': sigma_0_topo, 'system_name': 'both_sc',
                         'size': plasmid_size, 'structure': 'circular',
                         'topoI_concentration': topoI_concentration, 'gyrase_concentration': gyrase_concentration,
                         'topoI': True, 'gyrase': True,
                         'topoI_bm': topoI_bm, 'topoI_em': topoI_em, 'topoI_unbm': topoI_unbm,
                         'gyrase_bm': gyrase_bm, 'gyrase_em': gyrase_em, 'gyrase_unbm': gyrase_unbm
                           }

    global_dict_both_rx = {'frames': frames, 'dt': dt, 'superhelicity': sigma_0_gyrase, 'system_name': 'both_rx',
                         'size': plasmid_size, 'structure': 'circular',
                         'topoI_concentration': topoI_concentration, 'gyrase_concentration': gyrase_concentration,
                         'topoI': True, 'gyrase': True,
                         'topoI_bm': topoI_bm, 'topoI_em': topoI_em, 'topoI_unbm': topoI_unbm,
                         'gyrase_bm': gyrase_bm, 'gyrase_em': gyrase_em, 'gyrase_unbm': gyrase_unbm
                        }

    # Create lists of conditions for each system
    # ------------------------------------------
    # Global dictionaries
    global_dict_list = [global_dict_topoI, global_dict_gyrase, global_dict_both_sc, global_dict_both_rx]

    # Arrays with global superhelical densities
    list_sigmas = [topoI_sigma, gyrase_sigma, both_sigma_sc, both_sigma_rx]

    # Finally, run objective function. run_objective_function will process our conditions
    # ----------------------------------------------------
    objective, exp1_output = pm.run_objective_function_exp1(global_dict_list=global_dict_list,
                                                             exp_superhelicals=list_sigmas,
                                                             n_simulations=n_simulations,
                                                             parallelization=parallel)
    # -----------------------------------------------------------------------------------------------
    # Experiment 2!
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

    topoI_site = Site(site_type='topoI', name='topoI_s', start=1000, end=0,
                 binding_model=copy.deepcopy(topoI_bm))

    # Global dictionaries - One for each system
    # ------------------------------------------
    global_dict_gyrase2 = {'frames': frames2, 'dt': dt, 'superhelicity': sigma0_RX, 'system_name': 'gyrase_sc',
                         'size': single_molecule_size, 'structure': 'linear',
                         'gyrase_concentration': gyrase_concentration,
                         'topoI': False, 'gyrase': True,
                         'gyrase_bm': None, 'gyrase_em': gyrase_em, 'gyrase_unbm': gyrase_unbm,
                         'sites': [gyrase_site]
                         }

    global_dict_topoI2 = {'frames': frames2, 'dt': dt, 'superhelicity': sigma0_NSC, 'system_name': 'topoI_rx',
                         'size': plasmid_size, 'structure': 'linear',
                         'topoI_concentration': topoI_concentration,
                         'topoI': True, 'gyrase': False, # If we want to include topoI and gyrase
                         'topoI_bm': None, 'topoI_em': topoI_em, 'topoI_unbm': topoI_unbm,
                         'sites': [topoI_site]
                         }

    # Run experiments
    # ----------------------------------------------------
    global_dict_list2 = [global_dict_gyrase2, global_dict_topoI2]

    exp2_output = pm.run_experiment_2(global_dict_list=global_dict_list2, n_simulations=n_simulations,
                                      parallelization=parallel)

    # Arrange outputs
    output_dict = {'loss': objective, 'status': STATUS_OK,  # hyperopt results
                   # Specific results for processing
                   'experiment_1_results': exp1_output,  # Kinetic experiment
                   'experiment_2_results': exp2_output}
    return output_dict

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Pre-Process
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Here, we build the reference curves from a kinetic model from Fenfei's paper:
#   Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching
# We use this to quantify the error of each test.

# -----------------------------------------------------
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
# Note that this sigmaf is not at the end of the simulation, but the sigma at which there is 0 Relaxed DNA.
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
topoI_sigma = sigma

# -----------------------------------------------------
# Build experimental curve for Gyrase
# -----------------------------------------------------
# Kinetics: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Product = Supercoiled DNA
# Substrate = Relaxed DNA; which initially is the same as the DNA concentration

# Integrate MM kinetics
# ------------------------------------------
# Initially, there's no supercoiled DNA, and all of the relaxed DNA concentration corresponds
# to the plasmid concentration.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_gyrase(vmax=v_max_gyrase, KM=K_M_gyrase,
                                                       Supercoiled_0=0.0, Relaxed_0=DNA_concentration,
                                                       frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
gyrase_sigma = sigma

# -----------------------------------------------------
# Build experimental curve for system with both Topo I and Gyrase (from Sc state)
# -----------------------------------------------------
# Kinetics Gyrase: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Kinetics Topoisomerase: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI

# Integrate MM kinetics
# ------------------------------------------
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_both_T_G(vmax_topoI=v_max_topoI, vmax_gyrase=v_max_gyrase,
                                                         KM_topoI=K_M_topoI, KM_gyrase=K_M_gyrase,
                                                         Supercoiled_0=DNA_concentration, Relaxed_0=0.0,
                                                         frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
both_sigma_sc = sigma

# -----------------------------------------------------
# Build experimental curve for system with both Topo I and Gyrase (from Rx state)
# -----------------------------------------------------
# Integrate MM kinetics
# ------------------------------------------
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_both_T_G(vmax_topoI=v_max_topoI, vmax_gyrase=v_max_gyrase,
                                                         KM_topoI=K_M_topoI, KM_gyrase=K_M_gyrase,
                                                         Supercoiled_0=0.0, Relaxed_0=DNA_concentration,
                                                         frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
both_sigma_rx = sigma

#----------------------------------------------------------------------------------------------------------------------
# Random Search algorithm!
#----------------------------------------------------------------------------------------------------------------------
trials = Trials() # We need this for hyperopt!

# Here, you define the parameter space for the random search. These are basically the limits and parameter sampling.
# For simplicity, we use a uniform distribution for sampling each parameter, and hopefully we explore the "complete"
# parameter space (which is multidimensional).
space = {

    # Topo I params
    # ----------------------------------
    # Binding
    'k_on_topoI': hp.uniform('k_on_topoI', 0.000042, 0.00042),
    'width_topoI': hp.uniform('width_topoI', 0.00125, 0.0125),
    'threshold_topoI': hp.uniform('threshold_topoI', -0.07, 0.03),
    # Effect
    'k_cat_topoI': hp.uniform('k_cat_topoI', 0.0, 20.0),
    # Unbinding
    'k_off_topoI': hp.uniform('k_off_topoI', 0.1, 1.0),

    # Gyrase params
    # ----------------------------------
    # Binding
    'k_on_gyrase': hp.uniform('k_on_gyrase', 0.000024, 0.00024),
    'width_gyrase': hp.uniform('width_gyrase', 0.00125, 0.0125),
    'threshold_gyrase': hp.uniform('threshold_gyrase', -0.06, 0.04),
    # Effect
    'k_cat_gyrase': hp.uniform('k_cat_gyrase', 0.0, 20.0),
    # We don't have a sense of these rates...
    'k_wrap_gyrase': hp.uniform('k_wrap_gyrase', 0.01, 1.0),
    'k_unwrap_gyrase': hp.uniform('k_unwrap_gyrase', 0.01, 1.0),
    'k_go_gyrase': hp.uniform('k_go_gyrase', 0.01, 1.0),
    'k_dwell_gyrase': hp.uniform('k_dwell_gyrase', 0.01, 1.0),
    # 'sigma0_gyrase': hp.uniform('sigma0_gyrase', -.15, -0.07),

    # Unbinding
    'k_off_gyrase': hp.uniform('k_off_gyrase', 0.1, 1.0)
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
for name in ['topoI', 'gyrase']:
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


