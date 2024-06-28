from TORCphysics import Circuit
import numpy as np
import multiprocessing
import pandas as pd
from TORCphysics import parallelization_tools as pt
from TORCphysics import topo_calibration_tools as tct


#  Create template circuit and sites. Change these conditions per system conditions. - You can create temporal
#       circuit and sites files per system. Then create new ones and overwrite.
#       Run in parallel each of the systems.
#       For each run, collect only data of interest (susceptibilities), or you could save as binary to study dynamics.
#       Quantities worth of getting would be production rates (unbinding), global supercoiling, supercoiling at promoter,
#       positions? Could be the case. Maybe KDEs like you did in RNAPTracking (so you don't accomulate too much data).

# **********************************************************************************************************************
# Inputs/Initial conditions
# **********************************************************************************************************************

# Simulation conditions
# --------------------------------------------------------------
#binding_model='MaxMinPromoterBinding'
binding_model='MaxMinPromoterBinding_cutoff'
outputf='production_rates'
dt = 0.25
initial_time = 0
final_time = 1000#3600 #9000 ~2.5hrs
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

# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 12 # 64

# Circuit conditions
# --------------------------------------------------------------
output_prefix = 'single_gene'  #Although we won't use it
series = True
continuation = False
enzymes_filename = None
environment_filename = 'environment.csv'
# We will define sites and circuit csv as we advance

# These two filenames will be created and updated for every single cases
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'


# Junier experimental data - we only need distances for now.
# --------------------------------------------------------------
weak_exp = '../junier_data/weak.csv'
medium_exp = '../junier_data/medium.csv'
strong_exp = '../junier_data/strong.csv'

# Promoter responses
# --------------------------------------------------------------
weak_response = '../promoter_responses/weak.csv'
medium_response = '../promoter_responses/medium.csv'
strong_response = '../promoter_responses/strong.csv'


# **********************************************************************************************************************
# Functions
# **********************************************************************************************************************
# These functions create csv files and write them so we can load them later for building Circuits

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
# Load data and prepare systems
# **********************************************************************************************************************

# Load exp
weak_exp = pd.read_csv(weak_exp)
medium_exp = pd.read_csv(medium_exp)
strong_exp = pd.read_csv(strong_exp)

# Extract distances
weak_distances = list(weak_exp['distance'])
medium_distances = list(medium_exp['distance'])
strong_distances = list(strong_exp['distance'])

# **********************************************************************************************************************
# Processs (In parallel)
# **********************************************************************************************************************

# Create a multiprocessing pool
pool = multiprocessing.Pool()

# Let's do it for different cases

# Weak promoter
# --------------------------------------------------------------
circuit_name = 'weak'
production_rates = []
for i, upstream_distance in enumerate(weak_distances):

    # Build circuit csv
    # --------------------------------------------------------------
    circuit_size = upstream_distance + gene_length + downstream_distance
    make_linear_circuit_csv(circuit_filename, circuit_size, sigma0, circuit_name)

    # Site csv
    # --------------------------------------------------------------
    start = upstream_distance
    end = upstream_distance + gene_length
    make_gene_site_csv(sites_filename, 'gene', 'reporter', start, end, 1,
                       binding_model, weak_response)


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
prod_rates['distance'] = weak_distances
prod_rates = pd.DataFrame.from_dict(prod_rates)
prod_rates.to_csv(outputf+'-'+circuit_name+'.csv', index=False)

# Medium promoter
# --------------------------------------------------------------
circuit_name = 'medium'
production_rates = []
for i, upstream_distance in enumerate(medium_distances):

    # Build circuit csv
    # --------------------------------------------------------------
    circuit_size = upstream_distance + gene_length + downstream_distance
    make_linear_circuit_csv(circuit_filename, circuit_size, sigma0, circuit_name)

    # Site csv
    # --------------------------------------------------------------
    start = upstream_distance
    end = upstream_distance + gene_length
    make_gene_site_csv(sites_filename, 'gene', 'reporter', start, end, 1,
                       binding_model, medium_response)

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
prod_rates['distance'] = medium_distances
prod_rates = pd.DataFrame.from_dict(prod_rates)
prod_rates.to_csv(outputf+'-'+circuit_name+'.csv', index=False)

# Strong promoter
# --------------------------------------------------------------
circuit_name = 'strong'
production_rates = []
for i, upstream_distance in enumerate(strong_distances):

    # Build circuit csv
    # --------------------------------------------------------------
    circuit_size = upstream_distance + gene_length + downstream_distance
    make_linear_circuit_csv(circuit_filename, circuit_size, sigma0, circuit_name)

    # Site csv
    # --------------------------------------------------------------
    start = upstream_distance
    end = upstream_distance + gene_length
    make_gene_site_csv(sites_filename, 'gene', 'reporter', start, end, 1,
                       binding_model, strong_response)

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
prod_rates['distance'] = strong_distances
prod_rates = pd.DataFrame.from_dict(prod_rates)
prod_rates.to_csv(outputf+'-'+circuit_name+'.csv', index=False)
