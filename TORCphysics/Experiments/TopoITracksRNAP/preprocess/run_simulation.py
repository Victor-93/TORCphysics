import numpy as np
import multiprocessing
import pandas as pd
from TORCphysics import parallelization_tools as pt

# Be sure to copy the topoisomerases csv files here!
# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 8#8#64#120 #60#48 #24#8 #96 # 120

# Simulation conditions
# --------------------------------------------------------------
dt = 1.0
#dt = 0.5
#dt = 0.25
initial_time = 0
final_time = 3600 #1hr
#final_time = 2000#500
time = np.arange(initial_time, final_time + dt, dt)

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = None
enzymes_filename = None
#environment_filename = 'environment_avg_dt'+str(dt)+'.csv'
environment_filename = 'noRNAP_environment.csv'
output_prefix = 'noRNAP'
frames = len(time)
series = True
continuation = False

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
frame_values = list(range(frames+1))

# Calculate the corresponding time values using dtposition_output
time_values = [dt * frame for frame in frame_values]

# Create a DataFrame with 'frames' and 'time' columns
time_df = pd.DataFrame({'frames': frame_values, 'time': time_values})

# Names to filter
enzymes_names = ['topoI', 'gyrase']

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

    cout = 'position_' + name + '_dt'+str(dt)+'.txt'
    np.savetxt(cout, pos_out)

#    cout = 'position_' + name + '_df.csv'
#    position_output_df.to_csv(cout, index=False, sep=',')

pool.close()
