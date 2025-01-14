import numpy as np
import pandas as pd
import sys
# sys.path.append("/users/ph1vvb")
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import params
import pickle
from hyperopt import tpe, hp, fmin, Trials
import os

# **********************************************************************************************************************
# Description
# **********************************************************************************************************************
# Based on the experiments of assessing the impact of genome architecture... This script calibrates models in
# TORCPhysics that best describe the experimental results.
# It outputs values for the best set of params, and a trials object that contains general information for each.
# It varies gamma

# **********************************************************************************************************************
# Inputs/Initial conditions - At least the ones you need to change before each run
# **********************************************************************************************************************
#model_code='GB-Stages-trials-avgx1-01-'
model_code = 'genearch_V1-01'

gamma = 0.01

#promoter_cases = ['NPROMOTER']
#dt=NDT
#initial_time = 0
#final_time = 5400 #~1.5hrs
#k_weak = NKWEAK # This is the one inferred from the experiment
#tests = 1500  # number of tests for parametrization

#susceptibility_based = False
susceptibility_based = True

promoter_cases = ['weak']
#promoter_cases = ['medium']
#promoter_cases = ['strong']

dt=1.0
initial_time = 0
final_time = 500 #5400# ~1.3hrs
n_simulations =  4 # #50 #12
k_weak = 0.02#3#25 #25 # This is the one inferred from the experiment

# Junier experimental data - we only need distances for now.for i, promoter_case in enumerate(promoter_cases):
experimental_files = []
promoter_responses_files = []
susceptibility_files = []

for pcase in promoter_cases:
    experimental_files.append('../../junier_data/inferred-rate_kw' + str(k_weak) + '_'+ pcase + '.csv')
    susceptibility_files.append('../../junier_data/'+ pcase + '.csv')

    # Promoter responses
    promoter_responses_files.append('../../promoter_responses/' + pcase + '.csv')

if susceptibility_based:
    info_file = 'sus_'+model_code + promoter_cases[0] + '_gamma'+ str(gamma) + '_dt' + str(dt)
else:
    info_file = model_code + promoter_cases[0] + '_gamma'+ str(gamma) + '-kw' + str(k_weak) + '_dt' + str(dt)
file_out = info_file

# Parallelization conditions
# --------------------------------------------------------------
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
print('n_simulations', n_simulations)
print('hyperopt tests', tests)
print('Total number of simulations per test:', 12*n_simulations)


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
environment_filename = 'environment_avgx1_dt'+str(dt)+'.csv'

# Load the promoter response
# -----------------------------------
presponse = pd.read_csv(promoter_responses_files[0])

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
                'gamma':gamma}

# RANGES FOR RANDOM SEARCH
# -----------------------------------
# Reporter
kon_min = 0.01
kon_max = .5
superhelical_op_min=-0.1
superhelical_op_max=0.0
spread_min=0.005
spread_max=0.05

# RNAP - rates
k_closed_min = 0.01
k_closed_max = 0.5
k_open_min = 0.01
k_open_max = 0.5
k_ini_min = 0.01  #- Maybe this k_ini we really don't need it
k_ini_max = 0.5
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
    # Choose the curve to use for optimization
    if susceptibility_based:
        exp_curve = susceptibility_curves
    else:
        exp_curve = experimental_curves

    my_objective, output_dict = tct.gene_architecture_run_simple_odict(big_global_list, big_variation_list,
                                                                 exp_curve, parallel_info,
                                                                 additional_results=additional_results,
                                                                 susceptibility_based=susceptibility_based)

    # Return objective
    if calibrating:
        return my_objective
    else:
        return my_objective, output_dict


# **********************************************************************************************************************
# Prepare system
# **********************************************************************************************************************

experimental_curves = []
susceptibility_curves = []
distances = []
files_list = []
promoter_responses = []
for i, promoter_case in enumerate(promoter_cases):

    exp_df = pd.read_csv(experimental_files[i])
    sus_df = pd.read_csv(susceptibility_files[i])

    # Load exp
    experimental_curves.append(exp_df)
    susceptibility_curves.append(sus_df)

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
# Optimization
# -----------------------------------------------------
trials = Trials()
space = {
    # SITE
    'kon': hp.uniform('kon', kon_min, kon_max),
    'spread': hp.uniform('spread', spread_min, spread_max),
    'superhelical_op': hp.uniform('superhelical_op', superhelical_op_min, superhelical_op_max),

    # RNAP
    'k_closed': hp.uniform('k_closed', k_closed_min, k_closed_max),
    'k_open': hp.uniform('k_open', k_open_min, k_open_max),
    'k_ini': hp.uniform('k_ini', k_ini_min, k_ini_max),
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
    print("Doing model " + model_code + '; case ' + str(promoter_cases[0]) + '; k_weak=' + str(k_weak) + '; dt=' + str(
        dt))
    print('Doing parallelization process for:')
    print('n_simulations', n_simulations)
    print('hyperopt tests', tests)
    print('Total number of simulations per test:', 12 * n_simulations)

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

# Let's save the trials object
# --------------------------------------------------------------------------
with open(file_out+'-trials.pkl', 'wb') as file:
    pickle.dump(trials, file)

# Save model!
# --------------------------------------------------------------------------

# This one will have all the values of gene architecture (SITE + RNAP) - luckly they don't have same params names!
complete_df = pd.DataFrame(columns=['k_on', 'superhelical_op', 'spread',
                                    'k_closed', 'k_open', 'width', 'threshold',
                                    'k_ini', 'gamma', 'kappa', 'velocity', 'stall_torque',
                                    'k_off'])
complete_df['k_on'] = best_df['kon']
complete_df['superhelical_op'] = best_df['superhelical_op']
complete_df['spread'] = best_df['spread']
complete_df['k_closed'] = best_df['k_closed']
complete_df['k_open'] = best_df['k_open']
complete_df['width'] = reporter_oparams['width']
complete_df['threshold'] = reporter_oparams['threshold']
complete_df['k_ini'] = best_df['k_ini']
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
        'k_on': va['kon'], 'superhelical_op': va['superhelical_op'],
        'spread': va['spread'], 'k_closed': va['k_closed'], 'k_open': va['k_open'],
        'width': reporter_oparams['width'], 'threshold': reporter_oparams['threshold'],
        'k_ini': va['k_ini'], 'gamma': RNAP_oparams['gamma'], 'kappa': RNAP_oparams['kappa'],
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
