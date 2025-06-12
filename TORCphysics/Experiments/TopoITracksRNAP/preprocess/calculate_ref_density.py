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
dt = 1.0
#dt = 0.5
#dt = 0.25
initial_time = 0
final_time = 3600 #1hr
#final_time = 200 #2000#500 #1000 - doesn't matter much
time = np.arange(initial_time, final_time + dt, dt)
#file_out = 'reference_' #+ name + '_dt' + str(dt) + '.txt'
file_out = 'avg_reference_'  # This one is for the environment where we use averaged values

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = None
enzymes_filename = None
environment_filename = 'noRNAP_environment.csv'
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
    x = np.loadtxt('position_'+name+'_dt'+str(dt)+'.txt')
    x = x[~np.isnan(x)]  # Just in case to remove nans

    # Calculate kde
    kde_x, kde_y = tct.calculate_KDE(x, nbins, scaled=True)

    # Let's save the kde_x and kde_y as reference.
    kde = np.column_stack((kde_x, kde_y))

    # Save the combined array to a text file
    np.savetxt(file_out+name+'_dt'+str(dt)+'.txt', kde)


