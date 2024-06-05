from TORCphysics import Circuit
import numpy as np
import multiprocessing
import pandas as pd
from TORCphysics import parallelization_tools as pt

# TODO: Problem with transcribing from 85 to -1000 for antitet
# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 12#96 # 120

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = '../enzymes.csv'
environment_filename = '../environment.csv'
output_prefix = 'WT_system'
frames = 5000 #50000
series = True
continuation = False
dt = 0.25

# Prepare parallelization
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

# Run parallelization and collect relevant data
# --------------------------------------------------------------
# Create a multiprocessing pool
pool = multiprocessing.Pool()
pool_results = pool.map(pt.single_simulation_return_dfs, Items)

# Process simulations outputs
# --------------------------------------------------------------
# Create a range of values for the 'frames' column
frame_values = list(range(frames+1))

# Calculate the corresponding time values using dt
time_values = [dt * frame for frame in frame_values]

# Create a DataFrame with 'frames' and 'time' columns
time_df = pd.DataFrame({'frames': frame_values, 'time': time_values})

# Names to filter
site_names = ['circuit', 'tetA', "mKalama1", 'mRaspberry', 'antitet', 'tetA']
enzyme_names = ['RNAP', 'topoI', 'gyrase']

# And filter site results starting by name
for name in site_names:

    print(name)
    # Prepare output dataframe
    superhelical_output_df = time_df.copy()
    binding_output_df = time_df.copy()
    number_enzymes_df = time_df.copy()

    # Extract info from dataframes
    for n, out_dict in enumerate(pool_results):
        enzymes_df = out_dict['enzymes_df']
        sites_df = out_dict['sites_df']
        environmental_df = out_dict['environmental_df']

        # Filter superhelical density
        if name == 'circuit':
            mask = sites_df['type'] == name
        else:
            mask = sites_df['name'] == name

        # Superhelical
        # ------------------------------------------------------------------------------
        superhelical = sites_df[mask]['superhelical'].to_numpy()
        simulation_df = pd.DataFrame({'simulation_'+str(n): superhelical})
        superhelical_output_df = pd.concat([superhelical_output_df, simulation_df], axis=1).reset_index(drop=True)

        # Binding Event
        # ------------------------------------------------------------------------------
        binding_event = sites_df[mask]['binding'].to_numpy()
        simulation_df = pd.DataFrame({'simulation_'+str(n): binding_event})
        binding_output_df = pd.concat([binding_output_df, simulation_df], axis=1).reset_index(drop=True)

        # #Enzymes
        # ------------------------------------------------------------------------------
        number_enzymes = sites_df[mask]['#enzymes'].to_numpy()
        simulation_df = pd.DataFrame({'simulation_'+str(n): number_enzymes})
        number_enzymes_df = pd.concat([number_enzymes_df, simulation_df], axis=1).reset_index(drop=True)


    cout = 'superhelical_' + name + '_df.csv'
    superhelical_output_df.to_csv(cout, index=False, sep=',')

    cout = 'binding_' + name + '_df.csv'
    binding_output_df.to_csv(cout, index=False, sep=',')

    cout = 'N_enzymes_' + name + '_df.csv'
    number_enzymes_df.to_csv(cout, index=False, sep=',')


# And now for the RNAPs
for name in enzyme_names:

    print(name)
    # Prepare output dataframe
    all_positions = np.empty([1])

    # Extract info from dataframes
    for n, out_dict in enumerate(pool_results):
        enzymes_df = out_dict['enzymes_df']
        sites_df = out_dict['sites_df']
        environmental_df = out_dict['environmental_df']

        # Filter
        mask = enzymes_df['name'] == name

        # Position
        # ------------------------------------------------------------------------------
        position = enzymes_df[mask]['position'].to_numpy()

        all_positions = np.concatenate([position, all_positions])


    cout = 'positions_' + name + '.txt'
    np.savetxt(cout, all_positions)

pool.close()
