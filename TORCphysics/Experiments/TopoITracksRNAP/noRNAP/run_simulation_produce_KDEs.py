import numpy as np
import multiprocessing
import pandas as pd
from TORCphysics import parallelization_tools as pt
from TORCphysics import Circuit
from TORCphysics import topo_calibration_tools as tct
import matplotlib.pyplot as plt

# WAIT: At the moment, it doesn't do much
# Be sure to copy the topoisomerases csv files here!
# --------------------------------------------------------------
n_simulations = 4#8#8#64#120 #60#48 #24#8 #96 # 120
sets = 4

# Simulation conditions
# --------------------------------------------------------------
dt = 1.0 #0.25
initial_time = 0
final_time = 500
time = np.arange(initial_time, final_time + dt, dt)

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = None
enzymes_filename = None
environment_filename = 'environment_dt'+str(dt)+'.csv'
output_prefix = 'noRNAP'
frames = len(time)
#frames = 1000#5000 #50000
series = True
continuation = False
#dt = 0.25 #1#0.5 #  0.25

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

names = ['topoI', 'gyrase']
colors_dict = {'topoI': 'red', 'gyrase': 'cyan'}
kwargs = {'linewidth': 2, 'ls': '-'}

# Prepare parallelization
# --------------------------------------------------------------
global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
               'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
               'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
               'frames': frames, 'dt': dt, 'n_simulations': n_simulations}

Items = []
s = 0
for n in range(sets):
    item_small = []
    for simulation_number in range(n_simulations):
        g_dict = dict(global_dict)
        g_dict['n_simulations'] = s
        item_small.append(g_dict)
        s=s+1
    Items.append(item_small)

# Prepare processing info
# --------------------------------------------------------------
# Let's load the circuit so we have some info
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

frame_values = list(range(frames+1))

# Calculate the corresponding time values using dtposition_output
time_values = [dt * frame for frame in frame_values]

# Create a DataFrame with 'frames' and 'time' columns
time_df = pd.DataFrame({'frames': frame_values, 'time': time_values})

# Names to filter
enzymes_names = ['topoI', 'gyrase']

KDE_dict = {name: [] for name in enzymes_names}


# Run parallelization and produce KDEs
# --------------------------------------------------------------
# Create a multiprocessing pool
pool = multiprocessing.Pool(n_simulations)

for n in range(sets):
    pool_results = pool.map(pt.single_simulation_return_dfs, Items[n])

    # Process simulations outputs
    # --------------------------------------------------------------

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

            # Append positions
            # ------------------------------------------------------------------------------
            position = enzymes_df[mask]['position'].to_numpy()
            pos_out = np.concatenate([position, pos_out])


        # KDEs stuff
        # ------------------------------------------------------------------------------
        # number of bins for histogram - which we will not plot
        nbins = tct.calculate_number_nbins(my_circuit, name)

        # Calculate KDE
        kde = tct.calculate_KDE(data=pos_out, nbins=nbins, scaled=True)

        KDE_dict[name].append(kde)

# Calculate average, plot, and save
# --------------------------------------------------------------
fig, axs = plt.subplots(2, figsize=(width, 2*height), tight_layout=True)
for p, name in enumerate(names):

    # Processing stuf
    # -----------------------------------------------

    kde_x = KDE_dict[name][0][0]
    # Extract the y values (the second element in each tuple)
    y_values = [t[1] for t in KDE_dict[name][:]]

    # Calculate the mean and standard deviation
    kde_y = np.mean(y_values, axis=0)
    kde_sy = np.std(y_values, axis=0)

    # Plotting stuff
    # -----------------------------------------------
    ax = axs[p]

    ax.plot(kde_x[:], kde_y[:], '-', lw=lw*.5, color='green', label='kde')
    #ax.fill_between(kde_x[:], kde_y[:]-kde_sy[:], kde_y[:]+kde_sy[:], '-', lw=lw*.5, color='green', label='kde', alpha=0.5)

    # Labels
    # ------------------------
    ax.set_ylabel('Density', fontsize=font_size)
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    ax.set_title(name, fontsize=title_size)

plt.show()
pool.close()
