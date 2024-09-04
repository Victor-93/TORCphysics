import numpy as np
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import params
import pickle
from hyperopt import tpe, hp, fmin, Trials
import os
import sys

# **********************************************************************************************************************
# Description
# **********************************************************************************************************************
# This time we want to calibrate using the inferred rates rather than the susceptibility (it is too volatile), and
# using the Spacer binding model + the RNAPStages effect model. We will optimize to obtain the rates that optimize
# the inferred rates curve (previously calculated by a pre-process). In this sense, we can do independent calibrations
# between promoters.

# **********************************************************************************************************************
# Inputs/Initial conditions - At least the ones you need to change before each run
# **********************************************************************************************************************

#promoter_cases = ['weak']
#promoter_cases = ['medium']
promoter_cases = ['strong']

# Junier experimental data - we only need distances for now.for i, promoter_case in enumerate(promoter_cases):
experimental_files = []
promoter_responses_files = []
for pcase in promoter_cases:
    experimental_files.append('../../junier_data/inferred-rate_' + pcase + '.csv')

    # Promoter responses
    promoter_responses_files.append('../../promoter_responses/' + pcase + '.csv')

info_file = 'Stages-' + promoter_cases[0] + '_dt1'
file_out = info_file

dt = 1.0 #0.25
initial_time = 0
final_time = 2000#5000  #30000 ~8.3hrs

# Parallelization conditions
# --------------------------------------------------------------
# So, we have 12 distances, we can do  6 sets, each set processes 2 distances? Note that some distances take more
# time to simulate, but it's ok.
# So, n_sets = 6, n_subsets = 1 or 2, because if we have 64 workers and n_innerworkers = 9, then we'll have 9-18 simulations
# per distance, which would be enough, right?

n_workers = 12  #64  # Total number of workers (cpus)
n_sets = 4  #4  # Number of outer sets
n_inner_workers = n_workers // (n_sets + 1)  # Number of workers per inner pool
n_subsets = 2 #1  #2  # Number of simulations per inner pool
# +1 because one worker is spent in creating the outer pool
tests = 10 # 100  # number of tests for parametrization

# Basically, you make 'tests' number of tests. Then, according to the number of distances (12), you create a number
# 'n_sets' of groups (or sets) to distribute the work. If you make 6 sets, then each set runs 2 different systems (distances)
# For each set, you have n_inner_workers assign to run simulations. And each inner worker will ron 'n_subsets' number of simulations.
# So is like:
# Outer pool creates: n_sets
# Inner pool distributed to: n_inner_workers
# Each n_inner_worker performs: n_subsets simulations
# Total number of simulations per test is: n_sets * n_subsets * n_inner_workers

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
environment_filename = 'environment_small.csv'

# Load the promoter response
# -----------------------------------
presponse = pd.read_csv(promoter_responses_files[0])

# Models to calibrate
# -----------------------------------
# Site
reporter_name = 'reporter'
reporter_type = 'site'
reporter_binding_model_name = 'SpacerBinding'
reporter_oparams = {'superhelical_op': -0.06, 'spacer': 17}  # Parameters for the spacer length.
# I believe that the three promoters (weak, medium & strong) have a spacer length of 17 nt.

# RNAP
k_ini = 0.2 # TODO: Modify this once you find out in the RNAPTrack parametrization
gamma = 0.8 #  And this one as well
RNAP_name = 'RNAP'
RNAP_type = 'environmental'
RNAP_effect_model_name = 'RNAPStagesStall'
RNAP_unbinding_model_name = 'RNAPStagesSimpleUnbinding'
RNAP_oparams = {'width': presponse['width'], 'threshold': presponse['threshold'], # This is related to the site...
                'velocity': params.v0, 'kappa': params.RNAP_kappa, 'stall_torque': params.stall_torque,
                'gamma':gamma, 'k_ini':k_ini}

# RANGES FOR RANDOM SEARCH
# -----------------------------------
# Reporter
spacer_kon_min = 0.005
spacer_kon_max = .5

# RNAP - rates
k_closed_min = 0.01
k_closed_max = 0.5
k_open_min = 0.01
k_open_max = 0.5
#k_ini_min = 0.05  - Maybe this k_ini we really don't need it
#k_ini_max = 0.5
k_off_min = 0.01
k_off_max = 0.5


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
                           'n_simulations': n_inner_workers,
                           'DNA_concentration': 0.0}
            global_list.append(global_dict)

        big_global_list.append(global_list)

        # Variation dictionaries
        # ------------------------------------------
        # Site

        binding_oparams = {'k_on': params['spacer_kon'], 'spacer': reporter_oparams['spacer'],
                           'superhelical_op': reporter_oparams['superhelical_op'], }

        reporter_variation = {'name': reporter_name, 'object_type': reporter_type,
                              'binding_model_name': reporter_binding_model_name, 'binding_oparams': binding_oparams}

        # RNAP
        effect_oparams = {'k_closed': params['k_closed'], 'k_open': params['k_open'],
                          'width': RNAP_oparams['width'], 'threshold': RNAP_oparams['threshold'],
                          'k_ini': RNAP_oparams['k_ini'],
                          'velocity': RNAP_oparams['velocity'], 'gamma': RNAP_oparams['gamma'],
                          'stall_torque': RNAP_oparams['stall_torque'], 'kappa': RNAP_oparams['kappa']}

        unbinding_oparams = {'k_off': params['k_off']}

        RNAP_variation = {'name': RNAP_name, 'object_type': RNAP_type,
                          'effect_model_name': RNAP_effect_model_name, 'effect_oparams': effect_oparams,
                          'unbinding_model_name': RNAP_unbinding_model_name, 'unbinding_oparams': unbinding_oparams}

        big_variation_list.append([reporter_variation, RNAP_variation])

    # Info needed for the parallelization
    parallel_info = {'n_workers': n_workers, 'n_sets': n_sets,
                     'n_subsets': n_subsets, 'n_inner_workers': n_inner_workers}

    # Finally, run objective function.
    # ------------------------------------------
    if calibrating:
        additional_results = False
    else:
        additional_results = True
    my_objective, output_dict = tct.gene_architecture_calibration_nsets_rates(big_global_list, big_variation_list,
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
        circuit_filename = circuit_name + '-circuit.csv'
        sites_filename = circuit_name + '-sites.csv'

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
# Optimization
# -----------------------------------------------------
trials = Trials()
space = {
    # SITE
    'spacer_kon': hp.uniform('spacer_kon', spacer_kon_min, spacer_kon_max),

    # RNAP
    'k_closed': hp.uniform('k_closed', k_closed_min, k_closed_max),
    'k_open': hp.uniform('k_open', k_open_min, k_open_max),

    'k_off': hp.uniform('k_off', k_off_min, k_off_max)

}

# Save the current standard output
original_stdout = sys.stdout

# Define the file where you want to save the output
output_file_path = info_file + '.info'

# Open the file in write mode
with open(output_file_path, 'w') as f:
    # Redirect the standard output to the file
    sys.stdout = f  ##

    # Your code that prints to the screen
    print("Hello, this is the info file for the calibration of gene architecture.")
    print('n_workers', n_workers)
    print('n_sets', n_sets)
    print('n_subsets', n_subsets)
    print('n_inner_workers', n_inner_workers)
    print('hyperopt tests', tests)
    print('Total number of simulations per test:', n_sets * n_subsets * n_inner_workers)
    print('Total number of actual workers:', n_sets * (1 + n_inner_workers))

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

best_df = pd.DataFrame.from_dict([best])

# Save model!
# --------------------------------------------------------------------------

# This one will have all the values of gene architecture (SITE + RNAP) - luckly they don't have same params names!
complete_df = pd.DataFrame(columns=['k_on', 'superhelical_op', 'spacer',
                                    'k_closed', 'k_open', 'width', 'threshold',
                                    'k_ini', 'gamma', 'kappa', 'velocity', 'stall_torque',
                                    'k_off'])
complete_df['k_on'] = best_df['spacer_kon']
complete_df['superhelical_op'] = reporter_oparams['superhelical_op']
complete_df['spacer'] = reporter_oparams['spacer']
complete_df['k_closed'] = best_df['k_closed']
complete_df['k_open'] = best_df['k_open']
complete_df['width'] = RNAP_oparams['width']
complete_df['threshold'] = RNAP_oparams['threshold']
complete_df['k_ini'] = RNAP_oparams['k_ini']
complete_df['gamma'] = RNAP_oparams['gamma']
complete_df['kappa'] = RNAP_oparams['kappa']
complete_df['velocity'] = RNAP_oparams['velocity']
complete_df['stall_torque'] = RNAP_oparams['stall_torque']
complete_df['k_off'] = best_df['k_off']
complete_df.to_csv(file_out + '.csv', index=False, sep=',')

# Let's save trials info (params and loses)
# --------------------------------------------------------------------------
params_df = pd.DataFrame(columns=['test', 'loss'])
for n in range(tests):
    tdi = trials.trials[n]  # dictionary with results for test n
    lo = trials.trials[n]['result']['loss']  # loss
    va = trials.trials[n]['misc']['vals']  #values
    # Add a new row using append method
    new_row = pd.DataFrame({
        'test': n, 'loss': lo,
        'k_on': va['spacer_kon'], 'superhelical_op': reporter_oparams['superhelical_op'],
        'spacer': reporter_oparams['spacer'], 'k_closed': va['k_closed'], 'k_open': va['k_open'],
        'width': RNAP_oparams['width'], 'threshold': RNAP_oparams['threshold'],
        'k_ini': RNAP_oparams['k_ini'], 'gamma': RNAP_oparams['gamma'], 'kappa': RNAP_oparams['kappa'],
        'velocity': RNAP_oparams['velocity'],
        'stall_torque': RNAP_oparams['stall_torque'], 'k_off': va['k_off']
    })

    params_df = pd.concat([params_df, new_row], ignore_index=True)

params_df.to_csv(file_out + '-values.csv', index=False, sep=',')

# Let's run the function once more with the best params to produce the data so we can then just plot it.
# --------------------------------------------------------------------------
objective, output = objective_function(params=best, calibrating=False)

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
