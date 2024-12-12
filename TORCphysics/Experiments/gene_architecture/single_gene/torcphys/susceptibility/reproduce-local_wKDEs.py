import numpy as np
import pandas as pd
# sys.path.append("/users/ph1vvb")
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import params
import pickle
import os
import sys

# **********************************************************************************************************************
# Description
# **********************************************************************************************************************
# After the calibration of GaussianBinding model + Stages, this program
# uses the outputted parametrisation to the re-run the cases, producing pkl files for analysis.

# **********************************************************************************************************************
# Inputs/Initial conditions - At least the ones you need to change before each run
# **********************************************************************************************************************

#promoter_cases = ['weak', 'medium', 'strong']
promoter_cases = ['weak']
#promoter_cases = ['medium']
#promoter_cases = ['strong']

dt=1.0
initial_time = 0
final_time = 5400#10000#5000  #30000 ~8.3hrs
n_simulations =  16#16 #50 #12 # Per system

k_weak = 0.02#3#25 #25 # This is the one inferred from the experiment

model_code='sus_GB-Stages-avgx2-02-'

# Junier experimental data - we only need distances for now.for i, promoter_case in enumerate(promoter_cases):
experimental_files = []
promoter_responses_files = []
susceptibility_files = []
for pcase in promoter_cases:
    experimental_files.append('../../junier_data/inferred-rate_kw' + str(k_weak) + '_'+ pcase + '.csv')
    susceptibility_files.append('../../junier_data/'+ pcase + '.csv')

    # Promoter responses
    promoter_responses_files.append('../../promoter_responses/' + pcase + '.csv')

param_file = model_code + promoter_cases[0] + '_dt' + str(dt)+'.csv' # This one has the parametrisation from calibration

info_file = 'reproduce-'+model_code + promoter_cases[0] + '_dt' + str(dt) + '-wKDEs'
file_out = info_file

# Parallelization conditions
# --------------------------------------------------------------
# So, we have 12 distances, we can do  6 sets, each set processes 2 distances? Note that some distances take more
# time to simulate, but it's ok.
# So, n_sets = 6, n_subsets = 1 or 2, because if we have 64 workers and n_innerworkers = 9, then we'll have 9-18 simulations
# per distance, which would be enough, right?

n_workers = 12#64#12  #64  # Total number of workers (cpus)
n_sets = 4#6#4  #4  # Number of outer sets
n_inner_workers = n_workers // (n_sets + 1)  # Number of workers per inner pool
n_subsets = 2#2 #1  #2  # Number of simulations per inner pool
# +1 because one worker is spent in creating the outer pool
tests = 5#10 #100  # number of tests for parametrization

# Basically, you make 'tests' number of tests. Then, according to the number of distances (12), you create a number
# 'n_sets' of groups (or sets) to distribute the work. If you make 6 sets, then each set runs 2 different systems (distances)
# For each set, you have n_inner_workers assign to run simulations. And each inner worker will ron 'n_subsets' number of simulations.
# So is like:
# Outer pool creates: n_sets
# Inner pool distributed to: n_inner_workers
# Each n_inner_worker performs: n_subsets simulations
# Total number of simulations per test is: n_sets * n_subsets * n_inner_workers

print("Doing model " + model_code + '; case ' +str(promoter_cases[0])+'; k_weak=' + str(k_weak) +'; dt=' +str(dt))
print('Doing parallelization process for:')
print('n_workers', n_workers)
print('n_sets', n_sets)
print('n_subsets', n_subsets)
print('n_inner_workers', n_inner_workers)
print('hyperopt tests', tests)
print('Total number of simulations per test:', n_sets * n_subsets * n_inner_workers)
print('Total number of actual workers:', n_sets * (1 + n_inner_workers))


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
#environment_filename = 'environment_avgx2_dt'+str(dt)+'_RNAPtracking_off.csv'
environment_filename = 'environment_avgx2_dt'+str(dt)+'.csv'
RNAP_filename = '../../avgx2_RNAP_dt'+str(dt)+'.csv'

# Load the promoter response
# -----------------------------------
presponse = pd.read_csv(promoter_responses_files[0])
RNAPresponse = pd.read_csv(RNAP_filename) # And RNAP

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
                'gamma':RNAPresponse['gamma'].iloc[0]}

# Let's read the csv file
# -----------------------------------
calibration_params = pd.read_csv(param_file)

# Convert the DataFrame to a dictionary (list of dictionaries)
calibration_dict = calibration_params.to_dict(orient='records')[0]  # Get the first (and only) dictionary

# Rename the key 'k_on' to 'kon'
if 'k_on' in calibration_dict:
    calibration_dict['kon'] = calibration_dict.pop('k_on')

# Now calibration_dict will have the key 'kon' instead of 'k_on'
print('This is the calibration params: ', calibration_dict)


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


# Optimization function
#-----------------------------------------------------------------------------------------------------------------------
# If calibrating=True, then only the objective function is returned (used by hyperopt). If it is false, then it returns
# the outputs
def objective_function(params, calibrating=True):
    # We need to prepare the inputs.
    # At the moment, we only have one system.

    #  Let's give the data necesary to the function, and it will do all on it's own.
    #  Here we need to just organize the inputs.
    #  We have 3, and 12
    #  If calibration=False, then return dfs and stuff like that? Yes!/

    big_global_list = []  # Each entry is a list of global dicts
    big_variation_list = []  # We have a variation for each promoter case (not distances).
    # We need globals, variations and references. References can be in experimental_curves
    # So we just need to build global and variations
    for i, promoter_case in enumerate(promoter_cases):

        files = files_list[i]
        dists = distances[i]

        global_list = []

        for j, upstream_distance in enumerate(dists):
            circuit_file = files['circuit'][j]
            site_file = files['sites'][j]
            # Global dictionaries
            # ------------------------------------------
            global_dict = {'circuit_filename': circuit_file, 'sites_filename': site_file,
                           'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                           'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                           'frames': frames, 'dt': dt,
                           'n_simulations': n_simulations,
                           'DNA_concentration': 0.0}
            global_list.append(global_dict)

        big_global_list.append(global_list)

        # Variation dictionaries
        # ------------------------------------------
        # Site
#    'k_closed': 0.01, 'k_open': 0.01, 'k_off': 0.01, 'k_ini': 0.3,
#    'width': presponse['width'].iloc[0], 'threshold': presponse['threshold'].iloc[0]

        binding_oparams = {
            # GaussianBinding
            'k_on': params['kon'], 'spread': params['spread'], 'superhelical_op': params['superhelical_op'],
            # Stages
            'k_open': params['k_open'], 'width': reporter_oparams['width'], 'threshold': reporter_oparams['threshold'],
            'k_ini': params['k_ini'], 'k_closed': params['k_closed'], 'k_off': params['k_off']
        }

        reporter_variation = {'name': reporter_name, 'object_type': reporter_type,
                              'binding_model_name': reporter_binding_model_name, 'binding_oparams': binding_oparams}

        # RNAP
        effect_oparams = {'velocity': RNAP_oparams['velocity'], 'gamma': RNAP_oparams['gamma'],
                          'stall_torque': RNAP_oparams['stall_torque'], 'kappa': RNAP_oparams['kappa']}

        unbinding_oparams = {}

        RNAP_variation = {'name': RNAP_name, 'object_type': RNAP_type,
                          'effect_model_name': RNAP_effect_model_name, 'effect_oparams': effect_oparams,
                          'unbinding_model_name': RNAP_unbinding_model_name, 'unbinding_oparams': unbinding_oparams}

        big_variation_list.append([reporter_variation, RNAP_variation])

    # Info needed for the parallelization
    parallel_info = {'n_simulations': n_simulations}

    # Finally, run objective function.
    # ------------------------------------------
    if calibrating:
        additional_results = False
    else:
        additional_results = True
    my_objective, output_dict = tct.gene_architecture_run_simple_wKDEs(big_global_list, big_variation_list,
                                                                              experimental_curves, parallel_info,
                                                                              additional_results=additional_results)

    # Return objective
    if calibrating:
        return my_objective
    else:
        return my_objective, output_dict


# **********************************************************************************************************************
# Prepare system
# **********************************************************************************************************************

experimental_curves = []
distances = []
files_list = []
promoter_responses = []
for i, promoter_case in enumerate(promoter_cases):

    exp_df = pd.read_csv(experimental_files[i])

    # Load exp
    experimental_curves.append(exp_df)

    # Extract distances
    dists = list(exp_df['distance'])
    distances.append(dists)

    flist = {'circuit': [], 'sites': []}
    for j, upstream_distance in enumerate(dists):
        # Produce csv files to load
        circuit_name = promoter_case + '_' + str(upstream_distance)

        # These two filenames will be created and updated for every single cases
        circuit_filename = circuit_name + model_code+ 'circuit_kw'+str(k_weak)+'_dt'+str(dt)+'.csv'
        sites_filename = circuit_name + model_code+'sites_kw'+str(k_weak)+'_dt'+str(dt)+'.csv'

        # Build circuit csv
        # --------------------------------------------------------------
        circuit_size = upstream_distance + gene_length + downstream_distance
        make_linear_circuit_csv(circuit_filename, circuit_size, sigma0, circuit_name)

        # Site csv
        # --------------------------------------------------------------
        start = upstream_distance
        end = upstream_distance + gene_length
        make_gene_site_csv(sites_filename, 'gene', 'reporter', start, end, 1,
                           reporter_binding_model_name, 'none')

        flist['circuit'].append(circuit_filename)
        flist['sites'].append(sites_filename)

    files_list.append(flist)

    presponse = pd.read_csv(promoter_responses_files[i])  #.to_dict()
    promoter_responses.append(presponse)


# **********************************************************************************************************************
# Let's run the function once more with the best params to produce the data so we can then just plot it.
# --------------------------------------------------------------------------
objective, output = objective_function(params=calibration_dict, calibrating=False)

# Save the dictionary to a file
with open(file_out + '.pkl', 'wb') as file:
    pickle.dump(output, file)

# **********************************************************************************************************************
# Let's erase files
# ---------------------------------------------------------------------------
for i, promoter_case in enumerate(promoter_cases):

    files = files_list[i]
    dists = distances[i]

    for j, upstream_distance in enumerate(dists):
        os.remove(files['circuit'][j])
        os.remove(files['sites'][j])
