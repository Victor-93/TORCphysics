import numpy as np
from TORCphysics import Circuit
import random
import sys
import pickle
from custom_gyrase import GyraseCyclesEffect, GyraseCyclesUnbinding

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This is an example for scripting, so we won't be saving data as csv files and instead put every data of interest
# into a dictionary and save it with pickle

# Let's try to test a plasmid of the same size as the one used in the Topoisomerase study
# with one binding site.

# Let's use a custom model for gyrase called GyraseCycles, and let's see how it performs visually.

# We want to run multiple simulations of three different scenarios:
# 1) Gyrase acting on relaxed DNA
# 2) Gyrase acting on positively supercoiled DNA
# 3) Topoisomerase I acting on negatively supercoiled DNA

# And we want to analyse how much twist they add throughout the simualtions

# ----------------------------------------------------------------------------------------------------------------------
# WHAT NEEDS DOING
# ----------------------------------------------------------------------------------------------------------------------
#TODO:
# 1. Let's variate parameters and see how the model performs.
# 2. The best way to do it might be with loops.

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0 #0.25
initial_time = 0
final_time = 10800
#final_time = 10800 // 10
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'twist-contribution'
n_sims = 4#100  # This is the total
n_sim_shown = 4

# For the simulation
circuit_filename = 'circuit_linear.csv' # This is like the one they used in the experiment
sites_gyrase = 'sites_gyrase.csv'  # This one has a gyrase binding site at the ~middle
sites_topoI = 'sites_topoI.csv'  # And this one a topoI binding site
enzymes_filename = None  # 'enzymes_test.csv'
environment_topoI = 'environment_topoI.csv'  # Only topoI
environment_gyrase = 'environment_gyrase.csv'  # Only gyrase

# Superhelical values (sigma) for each case.
sigma_0_topo = -0.11#0.075  # Approximately -20 supercoils according the paper
sigma_0_gyrase_relaxed = 0.0  # We suppose this one.
sigma_0_gyrase_positive = 0.11  # We suppose this one.

output_prefix = 'test0'
series = True
continuation = False

# Variation params
gyrase_oparams = {'k_on': 0.0002, 'k_off':0.25, 'k_wrap': 0.5, 'k_unwrap': 0.5, 'k_go':0.5, 'k_cat': 20., 'k_dwell':.25}

# ----------------------------------------------------------------------------------------------------------------------
# We assign the effect and unbinding model.
# For the Binding Model, we use the one we already had, the GyraseRecognition Model (this may change in the future).
def assign_GyraseCycles_models(circuit, **oparams):

    # Binding Model could be assigned in a similar way, but let's not do it yet.

    # Assign effect model
    gyrase_environment = [d for d in circuit.environmental_list if d.name == 'gyrase'][0]
    gyrase_effect_model = GyraseCyclesEffect(**oparams)
    gyrase_environment.effect_model = gyrase_effect_model
    gyrase_environment.effect_model_name = gyrase_effect_model.__class__.__name__
    gyrase_environment.effect_model_oparams = gyrase_effect_model.oparams
    gyrase_environment.effect_model_oparams = None

    # Unbinding Model
    gyrase_environment = [d for d in circuit.environmental_list if d.name == 'gyrase'][0]
    gyrase_unbinding_model = GyraseCyclesUnbinding(**oparams)
    gyrase_environment.unbinding_model = gyrase_unbinding_model
    gyrase_environment.unbinding_model_name = gyrase_unbinding_model.__class__.__name__
    gyrase_environment.unbinding_model_oparams = gyrase_unbinding_model.oparams
    gyrase_environment.unbinding_model_oparams = None

# ----------------------------------------------------------------------------------------------------------------------
# Simulations
# ----------------------------------------------------------------------------------------------------------------------

# For gyrase Relaxed ------------------------------------
environment_filename = environment_gyrase
sigma0 = sigma_0_gyrase_relaxed
enzymes_df_list = []
sites_df_list = []
sites_filename = sites_gyrase

for n in range(n_sims):
    # Initialize circuit with the initial conditions
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt)

    # Adding custom model section!!!!!!!!
    # ---------------------------------------------------------------------------
    # Get gyrase environment from environmental_list
    assign_GyraseCycles_models(my_circuit, **gyrase_oparams)

    # ---------------------------------------------------------------------------
    my_circuit.name = my_circuit.name + '_' + str(n) # We can change the name of the circuit
    my_circuit.seed = my_circuit.seed + n + random.randrange(sys.maxsize)  # Just in case so each simulation has a different random number generator
    my_circuit.rng = np.random.default_rng(my_circuit.seed)

    # Change superhelical density manually to sigma0
    my_circuit.reset_circuit_superhelicity(sigma0)

    # Run simulations but storing dataframes on memory then adding them to the lists.
    enzymes_df, sites_df, environment_df = my_circuit.run_return_dfs() # Function run_return_dfs() returns the dataframes with the results of the simulation (it does not write CSV files).

    # Append dataframes to the lists.
    enzymes_df_list.append(enzymes_df)
    sites_df_list.append(sites_df)

gyrase_dict = {'sites_df': sites_df_list}  # We can save all the info we want

# For gyrase Positive ------------------------------------
environment_filename = environment_gyrase
sigma0 = sigma_0_gyrase_positive
enzymes_df_list = []
sites_df_list = []
sites_filename = sites_gyrase

for n in range(n_sims):
    # Initialize circuit with the initial conditions
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt)

    # Adding custom model section!!!!!!!!
    # ---------------------------------------------------------------------------
    # Get gyrase environment from environmental_list
    assign_GyraseCycles_models(my_circuit, **gyrase_oparams)

    my_circuit.name = my_circuit.name + '_' + str(n) # We can change the name of the circuit
    my_circuit.seed = my_circuit.seed + n + random.randrange(sys.maxsize)
    my_circuit.rng = np.random.default_rng(my_circuit.seed)
    # Change superhelical density manually to sigma0
    my_circuit.reset_circuit_superhelicity(sigma0)

    # Run simulations but storing dataframes on memory then adding them to the lists.
    enzymes_df, sites_df, environment_df = my_circuit.run_return_dfs() # Function run_return_dfs() returns the dataframes with the results of the simulation (it does not write CSV files).

    # Append dataframes to the lists.
    enzymes_df_list.append(enzymes_df)
    sites_df_list.append(sites_df)

gyrase_pos_dict = {'sites_df': sites_df_list}

# For topoI ----------------------------------
environment_filename = environment_topoI
sigma0 = sigma_0_topo
enzymes_df_list = []
sites_df_list = []
sites_filename = sites_topoI

for n in range(n_sims):
    # Initialize circuit with the initial conditions
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt)
    my_circuit.name = my_circuit.name + '_' + str(n) # We can change the name of the circuit
    my_circuit.seed = my_circuit.seed + n + random.randrange(sys.maxsize)
    my_circuit.rng = np.random.default_rng(my_circuit.seed)

    # Change superhelical density manually to sigma0
    my_circuit.reset_circuit_superhelicity(sigma0)

    # Run simulations but storing dataframes on memory then adding them to the lists.
    enzymes_df, sites_df, environment_df = my_circuit.run_return_dfs() # Function run_return_dfs() returns the dataframes with the results of the simulation (it does not write CSV files).

    # Append dataframes to the lists.
    enzymes_df_list.append(enzymes_df)
    sites_df_list.append(sites_df)

topoI_dict = {'sites_df': sites_df_list}

big_dict = {'gyrase_positive': gyrase_pos_dict, 'gyrase': gyrase_dict, 'topoI': topoI_dict}

# We can save python objects in binary with pickle.
with open('example_data_cycles.pkl', 'wb') as file:
    pickle.dump(big_dict, file)
