from TORCphysics import Circuit
import numpy as np
import multiprocessing
import pandas as pd
from TORCphysics import parallelization_tools as pt
from TORCphysics import topo_calibration_tools as tct
import pickle
from hyperopt import tpe, hp, fmin, Trials
import os
import sys

# TODO: We want to make a protocol where we calibrate some parameters:
#  Which will be the gamma injection of RNAP. And the promoter response k_max/k_min.
#  First, let's use the curve of the weak promoter
#  Then we can integrate all 3 promoters
#  Let's do it in batches or nsets to speed up the process. And keep the number of tests small?
#  Let's start by varying gamma, k_max and k_min, and let's see what kind of behaviiours we get

# We want to calibrate the promoter responses (rate and threshold)

# **********************************************************************************************************************
# Inputs/Initial conditions
# **********************************************************************************************************************

promoter_cases = ['weak']
response_multiplier = [1] #, 2, 4]  # Multiplier response, so we don't have to explore so many params

# Junier experimental data - we only need distances for now.for i, promoter_case in enumerate(promoter_cases):
for pcase in promoter_cases:
    experimental_files = ['../../junier_data/' + pcase + '.csv']

    # Promoter responses
    promoter_responses_files = ['../../promoter_responses/' + pcase + '.csv']

info_file = 'genearc-v1'
# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 16  #12 # 16

# So, we have 12 distances, we can do  6 sets, each set processes 2 distances? Note that some distances take more
# time to simulate, but it's ok.
# So, n_sets = 6, n_subsets = 1 or 2, because if we have 64 workers and n_innerworkers = 9, then we'll have 9-18 simulations
# per distance, which would be enough, right?

n_workers = 12  #64#12  # Total number of workers (cpus)
n_sets = 6  # Number of outer sets
n_subsets = 1  #5#7#2  # Number of simulations per set
n_inner_workers = n_workers // (n_sets + 1)  # Number of workers per inner pool
# +1 because one worker is spent in creating the outer pool
tests = 2  # number of tests for parametrization

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
outputf = 'production_rates'
dt = 2  #0.25
initial_time = 0
final_time = 500  #20000#3600 #9000 ~2.5hrs
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
# We will define sites and circuit csv as we advance

# Models to calibrate
# -----------------------------------
reporter_name = 'reporter'
reporter_type = 'site'
reporter_binding_model_name = 'MaxMinPromoterBinding'

RNAP_name = 'RNAP'
RNAP_type = 'enzyme'
RNAP_effect_model_name = 'RNAPStall'

v0 = 30.0 # bps
stall_torque = 12.0
kappa = 0.5

# RANGES FOR RANDOM SEARCH
# -----------------------------------
k_max_min = 0.01
k_max_max = 0.1
k_min_min = 0.001
k_min_max = 0.01
gamma_min = 0.0
gamma_max = 1.0


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

    # TODO: This time we'll need to build the global dicts?
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
        presponse = promoter_responses[i]


        global_list = []

        for j, upstream_distance in enumerate(dists):
            circuit_file = files[j]['circuit']
            site_file = files[j]['site']
            # Global dictionaries
            # ------------------------------------------
            global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
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
        mult = response_multiplier[i]
        binding_oparams = {'k_max':mult*params['k_max'], 'k_min':params['k_min'],
                           'a':float(presponse['a']), 'b':float(presponse['b']),
                           'threshold':float(presponse['threshold']), 'width':float(presponse['width'])}

        reporter_variation = {'name': reporter_name, 'object_type': reporter_type,
                           'binding_model_name': reporter_binding_model_name, 'binding_oparams': binding_oparams}

        # RNAP
        effect_oparams = {'velocity': v0, 'stall_torque': stall_torque, 'kappa': kappa,
                          'gamma': params['gamma']}

        RNAP_variation = {'name': RNAP_name, 'object_type': RNAP_type,
                              'effect_model_name': RNAP_effect_model_name, 'effect_oparams': effect_oparams}

        big_variation_list.append( [reporter_variation, RNAP_variation] )

    # Info needed for the parallelization
    parallel_info = {'n_workers': n_workers, 'n_sets': n_sets,
                     'n_subsets': n_subsets, 'n_inner_workers': n_inner_workers}


    # TODO: Create tct function that does that parallelization with nsets.
    #  It needs as inputs the
    #  * big_global_list, big_variation_list, the experimental curves, distances, parallel info.
    #  (maybe you can make the experimental curve just an array with values, not a df to make it
    #  more versatile).
    #  * You also need to make another function that processes the outputs, or maybe add a flag?
    #    so you can return a dict with all the info (binding,unbinding,supercoiling (global and site), and kdes??)

    # Finally, run objective function.
    # ------------------------------------------
    my_objective, output_dict = tct.single_case_RNAPTracking_calibration_nsets_2scheme(global_dict_list,
                                                                                       variations_list,
                                                                                       list_reference,
                                                                                       target_dict)
    #    return my_objective
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
                           'MaxMinPromoterBinding', promoter_responses_files[i])

        flist['circuit'].append(circuit_filename)
        flist['sites'].append(sites_filename)

    files_list.append(flist)

    presponse = pd.read_csv(promoter_responses_files[i]).to_dict()
    promoter_responses.append(presponse)

# **********************************************************************************************************************
# OPTIMIZATION# Let's erase files
# ---------------------------------------------------------------------------
for i, promoter_case in enumerate(promoter_cases):

    files = files_list[i]
    dists = distances[i]

    for j, upstream_distance in enumerate(dists):
        os.remove(files['circuit'][j])
        os.remove(files['sites'][j])


# **********************************************************************************************************************
# Optimization
# -----------------------------------------------------
trials = Trials()
space = {
    # SITE
    'k_max': hp.uniform('k_max', k_max_min, k_max_max),
    'k_min': hp.uniform('k_min', k_min_min, k_min_max),
    # RNAP
    'gamma': hp.uniform('gamma', gamma_min, gamma_max)
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

# Let's erase files
# ---------------------------------------------------------------------------
for i, promoter_case in enumerate(promoter_cases):

    files = files_list[i]
    dists = distances[i]

    for j, upstream_distance in enumerate(dists):
        os.remove(files['circuit'][j])
        os.remove(files['sites'][j])


# Create a multiprocessing pool
pool = multiprocessing.Pool()
# Let's erase files
# ---------------------------------------------------------------------------
for i, promoter_case in enumerate(promoter_cases):

    files = files_list[i]
    dists = distances[i]

    for j, upstream_distance in enumerate(dists):
        os.remove(files['circuit'][j])
        os.remove(files['sites'][j])




# Process
# --------------------------------------------------------------
circuit_name = promoter_case
production_rates = []
for i, upstream_distance in enumerate(distances):

    # Build circuit csv
    # --------------------------------------------------------------
    circuit_size = upstream_distance + gene_length + downstream_distance
    make_linear_circuit_csv(circuit_filename, circuit_size, sigma0, circuit_name)

    # Site csv
    # --------------------------------------------------------------
    start = upstream_distance
    end = upstream_distance + gene_length
    make_gene_site_csv(sites_filename, 'gene', 'reporter', start, end, 1,
                       'MaxMinPromoterBinding', promoter_response)

    # Prepare parallel inputs
    # --------------------------------------------------------------

    global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                   'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                   'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                   'frames': frames, 'dt': dt, 'n_simulations': n_simulations}

    Items = []
    for simulation_number in range(n_simulations):
        g_dict = dict(global_dict)
        g_dict['n_simulations'] = simulation_number
        Items.append(g_dict)

    # Run simulation in parallel returning dfs.
    # --------------------------------------------------------------
    pool_results = pool.map(pt.single_simulation_return_dfs, Items)

    # --------------------------------------------------------------
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt)
    production_rate_dict = tct.process_pools_calculate_production_rate(pool_results, ['reporter'], my_circuit)
    production_rates.append(production_rate_dict['reporter'])

prod_rates = {}
prod_rates['prod_rate'] = production_rates
prod_rates['distance'] = distances
prod_rates = pd.DataFrame.from_dict(prod_rates)
prod_rates.to_csv(outputf + '-' + circuit_name + '.csv', index=False)
