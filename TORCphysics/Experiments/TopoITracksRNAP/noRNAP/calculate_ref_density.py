import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from TORCphysics import Circuit
from TORCphysics import topo_calibration_tools as tct


# Description
# ---------------------------------------------------------
# I will process the simulations results and will calculate histograms,
# and will fit a kde to the histograms and will save it, so we can have a reference density for calculating
# fold enrichment and correlations.

# Inputs
# ---------------------------------------------------------

# Simulation conditions - Even though we won't run the simulation again
# --------------------------------------------------------------
dt = 1.0 #0.25
initial_time = 0
final_time = 1000 #2000#500 #1000
time = np.arange(initial_time, final_time + dt, dt)

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = None
enzymes_filename = None
environment_filename = 'environment_small.csv'
output_prefix = 'noRNAP'
frames = len(time)
series = True
continuation = False

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

# Let's load the circuit so we have some info
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

# Let's plot
# ---------------------------------------------------------
fig, axs = plt.subplots(2, figsize=(width, 2*height), tight_layout=True)
for p, name in enumerate(names):

    ax = axs[p]

    nbins = tct.calculate_number_nbins(my_circuit, name)

    # Load data
    x = np.loadtxt('position_'+name+'.txt')
    x = x[~np.isnan(x)]  # Just in case to remove nans

    # Plot histogram
    hist = sns.histplot(x, bins=nbins, ax=ax, color=colors_dict[name], label=name)

    # Calculate kde
    kde_x, kde_y = tct.calculate_KDE(x, nbins, scaled=True)

    # Calculate interpolation
    x_interpolated = tct.get_interpolated_x(1, my_circuit.size)
    y_interpolated = tct.get_interpolated_kde(kde_x, kde_y, x_interpolated)

    ax.plot(x_interpolated, y_interpolated, '-', color='purple', lw=lw, label='interpolated')
    ax.plot(kde_x, kde_y, '-', lw=lw*.5, color='green', label='kde')

    # Labels
    # ------------------------
    ax.set_ylabel('Density', fontsize=font_size)
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    ax.set_title(name, fontsize=title_size)


    # Let's save the kde_x and kde_y as reference.
    kde = np.column_stack((kde_x, kde_y))

    # Save the combined array to a text file
    np.savetxt('reference_'+name+'.txt', kde)

plt.show()

