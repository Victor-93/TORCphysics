from TORCphysics import Circuit
import numpy as np
import multiprocessing
import pandas as pd
from TORCphysics import parallelization_tools as pt
from TORCphysics import topo_calibration_tools as tct
import pickle

# We want to calibrate the promoter responses (rate and threshold)

# **********************************************************************************************************************
# Inputs/Initial conditions
# **********************************************************************************************************************

promoter_case = 'weak'

# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 16  #12 # 16

# Junier experimental data - we only need distances for now.
# --------------------------------------------------------------
experimental_file = '../junier_data/' + promoter_case + '.csv'

# Promoter responses
# --------------------------------------------------------------
promoter_response = '../promoter_responses/' + promoter_case + '.csv'

# Simulation conditions
# --------------------------------------------------------------
outputf = 'production_rates'
dt = 1.0 #0.25
initial_time = 0
final_time = 20000#20000  #3600 #9000 ~2.5hrs
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

# These two filenames will be created and updated for every single cases
circuit_filename = promoter_case + '-circuit.csv'
sites_filename = promoter_case + '-sites.csv'


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


def build_processing_dict(pool_list, gdict):
    processing_items = []
    for pool_results in pool_list:
        # Sort the processing info
        # --------------------------------------------------
        my_circuit = tct.load_circuit(gdict)

        # Get target site
        target_gene = [site for site in my_circuit.site_list if site.name == 'reporter'][0]
        RNAP_env = [environment for environment in my_circuit.environmental_list if environment.name == 'RNAP'][0]

        # Define x-axes
        x_system = tct.get_interpolated_x(1, my_circuit.size)
        x_gene = tct.get_interpolated_x(target_gene.start - RNAP_env.size, target_gene.end)

        processing_dict = {'circuit': my_circuit, 'x_gene': x_gene, 'x_system': x_system, 'site_name': 'reporter',
                           'additional_results': True}

        processing_items.append((pool_results, processing_dict))
    return processing_items


# **********************************************************************************************************************
# Load data and prepare systems
# **********************************************************************************************************************

# Load exp
experimental_curve = pd.read_csv(experimental_file)

# Extract distances
distances = list(experimental_curve['distance'])

# **********************************************************************************************************************
# Processs (In parallel)
# **********************************************************************************************************************

# Create a multiprocessing pool
pool = multiprocessing.Pool()

# Process
# --------------------------------------------------------------
circuit_name = promoter_case
production_rates = []
all_results = []  # We will pickle this
for i, upstream_distance in enumerate(distances):

    # Build circuit csv
    # --------------------------------------------------------------
    circuit_size = upstream_distance + gene_length + downstream_distance
    make_linear_circuit_csv(circuit_filename, circuit_size, sigma0, circuit_name)

    # Site csv
    # --------------------------------------------------------------
    start = upstream_distance + 30 # Because of the RNAP polymerase starts behind the gene
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

    # Other processing
    processing_items = build_processing_dict(pool_results, global_dict)
    output_pools = pool.map(tct.gene_architecture_process_pool, processing_items)
    all_results.append({'distance': upstream_distance, 'results': output_pools})


pool.close()
prod_rates = {}
prod_rates['prod_rate'] = production_rates
prod_rates['distance'] = distances
prod_rates = pd.DataFrame.from_dict(prod_rates)
prod_rates.to_csv(outputf + '-' + circuit_name + '.csv', index=False)

# Save the dictionary to a file
with open(outputf + '-' + circuit_name +'.pkl', 'wb') as file:
    pickle.dump(all_results, file)