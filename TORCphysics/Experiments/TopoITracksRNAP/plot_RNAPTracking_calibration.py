import numpy as np
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import matplotlib.pyplot as plt
import pickle

# TODO: Once you have your "true" calibration results, modify the script to load them.
# Description
# --------------------------------------------------------------
# Following the calibration process by executing calibrate_tracking_nsets_p2.py for both the stalling and uniform
# models, here I plot their results

# Simulation conditions
# --------------------------------------------------------------
dt = 1.0 #0.5
initial_time = 0
final_time = 200 # 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Inputs
# --------------------------------------------------------------

pickle_files = [
    'track-uniform/calibration_RNAPTracking_nsets_p2_small_dt'+str(dt)+'.pkl',
    #'track-uniform/calibration_RNAPTracking_nsets_p2_small_dt0.5.pkl',
    'track-StagesStall/calibration_RNAPTracking_nsets_p2_small_dt'+str(dt)+'.pkl'
]
output_prefix = 'test-RNAPStages-topoIRNAPtracking'
title = ['Uniform model', 'Stall model']

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = None
environment_filename = 'environment_dt'+str(dt)+'.csv'
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

alpha=0.3

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


# Let's plot
# ---------------------------------------------------------
ncases = len(pickle_files)
fig, axs = plt.subplots(ncases, figsize=(width, ncases*height), tight_layout=True)

for n, pickle_file in enumerate(pickle_files):


    # Plot line at 1.0
    #axs[n].plot([-200,my_circuit.size+200], [1,1], '--', 'black', lw=1.0)

    # Load the dictionary from the file
    # ---------------------------------------------------------
    with open(pickle_file, 'rb') as file:
        output = pickle.load(file)


    for p, name in enumerate(names):

        ax = axs[n]

        # Different process for RNAP or topos
        if name == 'RNAP':
            ax2 = ax.twinx()  # Different scale for RNAP
            ax = ax2
            x = x_gene
            y = output['results']['kde_gene'][name]['mean']
            ys = output['results']['kde_gene'][name]['std']
            # ax2 = ax.twinx()
            # ax2.plot(x_gene, output['results']['kde_gene'][name]['mean'], colors_dict[name], lw=lw, label=name)
        else:
            x = x_system
            y = output['results']['FE_curve'][name]['mean']
            ys = output['results']['FE_curve'][name]['std']
            # ax.plot(x_system, output['results']['FE_curve'][name]['mean'], color=colors_dict[name], lw=lw, label=name)
        ax.plot(x, y, color=colors_dict[name], lw=lw, label=name)
        ax.fill_between(x, y-ys, y+ys, color=colors_dict[name], alpha=alpha)
        print(pickle_file, name, np.mean(y))

    # Labels
    # ------------------------
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax2.set_ylabel(ylabel2, fontsize=font_size, color=colors_dict['RNAP'])
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    ax.set_xlim(0, my_circuit.size)
    ax.grid(True)
    #ax.set_ylim(.5,1.5)
    ax.set_xlim(0, my_circuit.size)
    ax.set_title(title[n], fontsize=font_size)
    # fig.suptitle(title[n], fontsize=title_size)
    if n == 0:
        ax.legend(loc='best', fontsize=font_size)


    # Add FE and correlation labels
    # -------------------------------------
    FE = output['FE'][0]
    CO = output['overall_correlation']
    OB = output['objective']
    RNAP_CO = output['results']['RNAP_correlation']

    # Define the text to display
    textstr = f'FE={FE:.2f}, RCO={RNAP_CO:.2f}, CO={CO:.2f}, OB={OB:.2e}'

    # Add the text box to the plot
    props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
    ax.text(0.05, 0.1, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    print('For ', title[n])
    print('FE',FE)
    print('CO',CO)
    print('objective', output['objective'])
    print('RNAP_correlation', output['results']['RNAP_correlation'])

plt.show()
#plt.savefig(output_prefix+'-FE.png')
#plt.savefig(output_prefix+'-FE.pdf')
