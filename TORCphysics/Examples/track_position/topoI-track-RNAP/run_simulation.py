import time

from TORCphysics import Circuit
import numpy as np
import multiprocessing
import pandas as pd
from TORCphysics import parallelization_tools as pt

# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 40 #96 # 120

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'topoIRNAPtrack'
frames = 5000 #50000
series = True
continuation = False
dt = 0.5 #  0.25

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
enzymes_names = ['RNAP', 'topoI', 'gyrase']

# And filter results starting by name
for name in enzymes_names:

    # Prepare output dataframe
    position_output_df = time_df.copy()
    pos_out = np.empty(1)

    # Extract info from dataframes
    for n, out_dict in enumerate(pool_results):
        enzymes_df = out_dict['enzymes_df']
        sites_df = out_dict['sites_df']
        environmental_df = out_dict['environmental_df']

        # Filter superhelical density
        mask = enzymes_df['name'] == name

        # Position
        # ------------------------------------------------------------------------------
        position = enzymes_df[mask]['position'].to_numpy()
        pos_out = np.concatenate([position, pos_out])
        # simulation_df = pd.DataFrame({'simulation_'+str(n): position})
        # position_output_df = pd.concat([position_output_df, simulation_df], axis=1).reset_index(drop=True)

    cout = 'position_' + name + '.txt'
    np.savetxt(cout, pos_out)

#    cout = 'position_' + name + '_df.csv'
#    position_output_df.to_csv(cout, index=False, sep=',')

pool.close()
