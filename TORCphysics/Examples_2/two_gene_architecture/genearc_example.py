import numpy as np
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import params
from TORCphysics import Circuit
import random
import pickle
import os
import sys

# **********************************************************************************************************************
# Description
# **********************************************************************************************************************
# Based on the experiments of assessing the impact of genome architecture...

# TODO: 1) Try adding a second gene and modifying the orientation. 2) Mixing promoters
#       3) Increasing the distances between the gene, and the barriers.
#

# **********************************************************************************************************************
# Inputs/Initial conditions - At least the ones you need to change before each run
# **********************************************************************************************************************
promoter_case = 'weak'
dt=1.0
initial_time = 0
final_time = 5400 #~1.5hrs

RNAP_gamma = 0.1

n_sims = 4
file_out = 'genearch_example'
# Simulation conditions
# --------------------------------------------------------------
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# The system is like:  UP______GENE____DOWN, where UP and DOWN are upstream and downstream barriers.
# The upstream barrier will always be located at 0 (linear system), the gene start site is located at
# x=upstream_distance and termination at upstream_distance+gene_length, and the total length of the region is
# upstream_distance + gene_length + downstream_distance
gene_length = 900
downstream_distance = 320

# Initial superhelical density
sigma0 = -0.046

# Circuit conditions
# --------------------------------------------------------------
output_prefix = 'single_gene'  #Although we won't use it
series = True
continuation = False
enzymes_filename = None
environment_filename = 'environment.csv'

# Load the promoter response
# -----------------------------------
promoter_responses_file = promoter_case + '.csv'
presponse = pd.read_csv(promoter_responses_file)

# Models to calibrate
# -----------------------------------
# Site
reporter_name = 'reporter'
reporter_type = 'site'
reporter_binding_model_name = 'GaussianBinding'
reporter_oparams = {
    # Parameters for the spacer length.
    'k_on': 0.01, 'superhelical_op': -0.06, 'spread': 0.05,
    # These are related to the site...
    'k_closed': 0.01, 'k_open': 0.01, 'k_off': 0.01, 'k_ini': 0.3,
    'width': presponse['width'].iloc[0], 'threshold': presponse['threshold'].iloc[0]
}

# RNAP
RNAP_name = 'RNAP'
RNAP_type = 'environmental'
RNAP_effect_model_name = 'RNAPStagesStallv2'
RNAP_unbinding_model_name = 'RNAPStagesSimpleUnbindingv2'
RNAP_oparams = {'velocity': params.v0, 'kappa': params.RNAP_kappa, 'stall_torque': params.stall_torque,
                'gamma': RNAP_gamma} # HERE WE SET THE VAL

# **********************************************************************************************************************
# Functions
# **********************************************************************************************************************
# These functions create csv files and write them so we can load them later for building Circuits
#-----------------------------------------------------------------------------------------------------------------------
# Makes a circuit.csv file g
def make_linear_circuit_csv(filename, nbp, sigma0, name):
    info = {'name': name, 'structure': 'linear', 'size': nbp, 'twist': 0.0, 'superhelical': sigma0, 'sequence': 'none'}
    circuit_df = pd.DataFrame([info])
    circuit_df.to_csv(filename, index=False)


# Makes a site csv containing one gene (could be expanded to more sites).
def make_gene_site_csv(filename, stype, name, start, end, k_on, bmodel, paramsfile):
    info = {
        'type': stype, 'name': name, 'start': start, 'end': end, 'k_on': k_on,
        'binding_model': bmodel, 'binding_oparams': paramsfile
    }
    site_df = pd.DataFrame([info])
    site_df.to_csv(filename, index=False)


# **********************************************************************************************************************
# Prepare system
# **********************************************************************************************************************

experimental_curves = []
susceptibility_curves = []
distances = []
files_list = []
promoter_responses = []

barrier_distance = 500

upstream_distance = barrier_distance

circuit_name = promoter_case + '_' + str(upstream_distance)

# These two filenames will be created and updated for every single cases
circuit_filename = circuit_name + '_circuit.csv'
sites_filename = circuit_name + '_sites.csv'

# Build circuit csv
# --------------------------------------------------------------
circuit_size = upstream_distance + gene_length + downstream_distance
make_linear_circuit_csv(circuit_filename, circuit_size, sigma0, circuit_name)

# Site csv
# --------------------------------------------------------------
start = upstream_distance
end = upstream_distance + gene_length
make_gene_site_csv(sites_filename, 'gene', promoter_case, start, end, 1,
                   reporter_binding_model_name, 'none')

# **********************************************************************************************************************
# RUN SIMULATION
# **********************************************************************************************************************
enzymes_df_list = []
sites_df_list = []
environment_df_list = []
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
    environment_df_list.append(environment_df)

results_dict = {'sites_df': sites_df_list, 'enzymes_df': enzymes_df_list, 'environment_df': environment_df_list}

# Save the dictionary to a file
with open(file_out+'.pkl', 'wb') as file:
    pickle.dump(results_dict, file)

# **********************************************************************************************************************
# OPTIONAL: Erase circuit files generated
# ---------------------------------------------------------------------------
os.remove(circuit_filename)
os.remove(sites_filename)
