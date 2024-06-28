from TORCphysics import Circuit
import numpy as np
import multiprocessing
import pandas as pd
from TORCphysics import parallelization_tools as pt
from TORCphysics import topo_calibration_tools as tct


# This time we want to get the DFs, of unbinding

# **********************************************************************************************************************
# Inputs/Initial conditions
# **********************************************************************************************************************

promoter_case = 'weak'

# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 4 # 16
idist =11# 0,4, 11


# Junier experimental data - we only need distances for now.
# --------------------------------------------------------------
experimental_file = '../../junier_data/' + promoter_case + '.csv'

# Promoter responses
# --------------------------------------------------------------
promoter_response = '../../promoter_responses/'+promoter_case+'.csv'

# Simulation conditions
# --------------------------------------------------------------
outputf='production_rates'
dt = 0.25
initial_time = 0
final_time = 20000#36000#9000#3600 #9000 ~2.5hrs
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
circuit_filename = promoter_case+'-circuit.csv'
sites_filename = promoter_case+'-sites.csv'

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
for i, upstream_distance in enumerate(distances):

    print(i, upstream_distance)

    if i != idist:
        continue
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

    # Process results
    # --------------------------------------------------------------
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt)

    # Get target site
    target_gene = [site for site in my_circuit.site_list if site.name == 'reporter'][0]
    RNAP_env = [environment for environment in my_circuit.environmental_list if environment.name == 'RNAP'][0]

    # Define x-axes
    x_system = tct.get_interpolated_x(1, my_circuit.size)
    x_gene = tct.get_interpolated_x(target_gene.start - RNAP_env.size, target_gene.end)

    # Create a range of values for the 'frames' column
    frame_values = list(range(my_circuit.frames + 1))

    # Calculate the corresponding time values using dt
    time_values = [my_circuit.dt * frame for frame in frame_values]

    # Create a DataFrame with 'frames' and 'time' columns
    time_df = pd.DataFrame({'frames': frame_values, 'time': time_values})

    # Let's do it for sites first!
    # ------------------------------------------------------------------------------
    # Prepare output dataframe
    superhelical_output_df = time_df.copy()
    binding_output_df = time_df.copy()
    unbinding_output_df = time_df.copy()
    nenzymes_output_df = time_df.copy()

    # And filter site results starting by name
    for name in ['reporter']:
        for n, out_dict in enumerate(pool_results):

            #Get output dfs
            sites_df = out_dict['sites_df']

            # Filter superhelical density
            if name == 'circuit':
                mask = sites_df['type'] == name
            else:
                mask = sites_df['name'] == name

            mask_circuit = sites_df['type'] == 'circuit'

            # Collect measurements
            # ------------------------------------------------------------------------------
            gsigma = sites_df[mask_circuit]['superhelical'].to_numpy()
            superhelical = sites_df[mask]['superhelical'].to_numpy()
            binding_event = sites_df[mask]['binding'].to_numpy()
            unbinding_event = sites_df[mask]['unbinding'].to_numpy()
            number_enzymes = sites_df[mask]['#enzymes'].to_numpy()

            # Make them DFs and sort by simulation
            gsigma = pd.DataFrame({'simulation_' + str(n): gsigma})
            superhelical = pd.DataFrame({'simulation_' + str(n): superhelical})
            binding_event = pd.DataFrame({'simulation_' + str(n): binding_event})
            unbinding_event = pd.DataFrame({'simulation_' + str(n): unbinding_event})
            number_enzymes = pd.DataFrame({'simulation_' + str(n): number_enzymes})

            if n == 0:
                gsigma_df = gsigma.copy()
                superhelical_df = superhelical.copy()
                binding_event_df = binding_event.copy()
                unbinding_event_df = unbinding_event.copy()
                number_enzymes_df = number_enzymes.copy()
            else:
                gsigma_df = pd.concat([gsigma_df, gsigma], axis=1).reset_index(drop=True)
                superhelical_df = pd.concat([superhelical_df, superhelical], axis=1).reset_index(drop=True)
                binding_event_df = pd.concat([binding_event_df, binding_event], axis=1).reset_index(drop=True)
                unbinding_event_df = pd.concat([unbinding_event_df, unbinding_event], axis=1).reset_index(drop=True)
                number_enzymes_df = pd.concat([number_enzymes_df, number_enzymes], axis=1).reset_index(drop=True)

    binding_event_df.to_csv(promoter_case+'_binding_distance-'+str(int(upstream_distance))+'.csv', index=False)
    gsigma_df.to_csv(promoter_case+'_gsigma_distance-'+str(int(upstream_distance))+'.csv', index=False)
    superhelical_df.to_csv(promoter_case+'_superhelical_distance-'+str(int(upstream_distance))+'.csv', index=False)
    unbinding_event_df.to_csv(promoter_case+'_unbinding_distance-'+str(int(upstream_distance))+'.csv', index=False)


    # For the enzymes.
    # ------------------------------------------------------------------------------
    enzymes_dict = {}
    for name in ['RNAP', 'topoI', 'gyrase']:

        all_positions = np.empty([1])

        # Extract info from dataframes
        for n, out_dict in enumerate(pool_results):
            enzymes_df = out_dict['enzymes_df']

            # Filter
            mask = enzymes_df['name'] == name

            # Position
            # ------------------------------------------------------------------------------
            position = enzymes_df[mask]['position'].to_numpy()

            all_positions = np.concatenate([position, all_positions])


        # number of bins for histogram - which we will not plot
        if name == 'RNAP':
            nbins = int(gene_length/50)
            size = abs(RNAP_env.site_list[0].start - RNAP_env.size - RNAP_env.site_list[0].end)
            nbins = int(size / RNAP_env.size)  # because genes are smaller
            nbins=90
            x = x_gene
        else:
            my_environmental = \
            [environmental for environmental in my_circuit.environmental_list if environmental.name == name][
                0]
            nbins = int(my_circuit.size / my_environmental.size)  # number of bins. - note that this only applies for topos

            #nbins = tct.calculate_number_nbins(my_circuit, name)
            x = x_system
            #nbins=20

        print(name, nbins)
        # Calculate KDE
        kde_x, kde_y = tct.calculate_KDE(data=all_positions, nbins=nbins, scaled=True)

        kde_y = tct.get_interpolated_kde(kde_x, kde_y, x)
        kde_x = x

        kde = np.column_stack((kde_x, kde_y))

        #kde = np.dstack((kde_x, kde_y))
        np.savetxt(promoter_case+'_'+name+'_KDE_distance-'+str(int(upstream_distance))+'.csv', kde)

    break

