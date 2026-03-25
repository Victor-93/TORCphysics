import numpy as np
from TORCphysics import Circuit
import random
import sys
import pickle

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# The idea of this script is to quickly produce a simulation with the standard gyrase paramterisation,
# so it can be compared with the calibration of the new cycles model.

# We want to run the following scenario: Gyrase acting on relaxed DNA

# And we want to analyse how much twist they add throughout the simualtions...

# Make sure to run the simulations with the same conditions.
# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0 #0.25
initial_time = 0
final_time = 10800
#final_time = 10800 // 2
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
environment_gyrase = 'environment_gyrase.csv'  # Only gyrase

# Superhelical values (sigma) for each case.
sigma_0_topo = -0.11#0.075  # Approximately -20 supercoils according the paper
sigma_0_gyrase_relaxed = 0.0  # We suppose this one.
sigma_0_gyrase_positive = 0.11  # We suppose this one.

output_prefix = 'test0'
series = True
continuation = False

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

# We can save python objects in binary with pickle.
with open('standard_gyrase-model.pkl', 'wb') as file:
    pickle.dump(gyrase_dict, file)
