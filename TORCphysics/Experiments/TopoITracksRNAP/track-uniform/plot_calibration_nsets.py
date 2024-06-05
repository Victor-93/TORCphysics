import numpy as np
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import matplotlib.pyplot as plt
import pickle

# Description
# --------------------------------------------------------------
# Following the calibration process by executing calibrate_tracking_nsets.py, this script loads the pickle output
# to produce the necessary plots for the analysis.

# Inputs
# --------------------------------------------------------------
pickle_file = 'calibration_RNAPTracking_nsets_p2.pkl'

# Simulation conditions
# --------------------------------------------------------------
dt = 0.25
initial_time = 0
final_time = 500 # o 2000
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = '../sites.csv'
enzymes_filename = None
environment_filename = 'environment_calibration.csv'
output_prefix = 'topoIRNAPtrack-uniform'
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

names = ['RNAP', 'topoI', 'gyrase']
colors_dict = {'RNAP': 'purple', 'topoI': 'red', 'gyrase': 'cyan'}
kwargs = {'linewidth': 2, 'ls': '-'}
ylabel = 'Fold Enrichment'
ylabel2 = 'RNAP Density'

# TARGETS FOR OPTIMIZATION - We need this even though we will not run the parallelization
# -----------------------------------
target_FE = 1.68  # Target fold-enrichment
target_CO = 1.0  # Target correlation between topo I and RNAP densities.
x_spacing = 10.0  # The spacing I want at the moment of doing interpolation.

# nbins is the number of bins to use when calculating the kde
target_dict = {'target_FE': target_FE, 'target_CO': target_CO, 'target_gene': 'reporter',
               'enzymes_names': ['RNAP', 'topoI', 'gyrase']}


# Let's load the circuit, so we can extract some information
# -----------------------------------
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

# Get target site
target_gene = [site for site in my_circuit.site_list if site.name == target_dict['target_gene']][0]
RNAP_env = [environment for environment in my_circuit.environmental_list if environment.name == 'RNAP'][0]

# Define x-axes
x_system = tct.get_interpolated_x(1, my_circuit.size)
x_gene = tct.get_interpolated_x(target_gene.start - RNAP_env.size, target_gene.end)


# Load the dictionary from the file
# ---------------------------------------------------------
with open(pickle_file, 'rb') as file:
    output = pickle.load(file)

# Let's plot
# ---------------------------------------------------------
fig, axs = plt.subplots(1, figsize=(width, 1*height), tight_layout=True)
for p, name in enumerate(names):

    ax = axs


    # Different process for RNAP or topos
    if name == 'RNAP':
        ax2 = ax.twinx()  # Different scale for RNAP
        ax2.plot(x_gene, output['results']['kde_gene'][name]['mean'], colors_dict[name], lw=lw, label=name)
    else:
        ax.plot(x_system, output['results']['FE_curve'][name]['mean'], color=colors_dict[name], lw=lw, label=name)

# Labels
# ------------------------
ax.set_ylabel(ylabel, fontsize=font_size)
ax2.set_ylabel(ylabel2, fontsize=font_size, color=colors_dict['RNAP'])
ax.set_xlabel(r'Position (bp)', fontsize=font_size)
ax.set_xlim(0, my_circuit.size)
ax.grid(True)
ax.legend(loc='best', fontsize=font_size)

# Add FE and correlation labels
# -------------------------------------
FE = output['FE'][0]
CO = output['overall_correlation']

# Define the text to display
textstr = f'FE={FE:.2f}, CO={CO:.2f}'

# Add the text box to the plot
props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
ax.text(0.05, 0.1, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
print('FE',FE)
print('CO',CO)
print('objective', output['objective'])
plt.show()
