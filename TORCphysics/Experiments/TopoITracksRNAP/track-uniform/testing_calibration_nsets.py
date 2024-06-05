import numpy as np
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import matplotlib.pyplot as plt

# Description
# --------------------------------------------------------------
# Very similar to testing_calinration.py, but now the porpuse of this script is to test the
# single_case_RNAPTracking_calibration_nsets function from topo_calibration_tools.py
# This new calibration protocal launches n simulations for a number of tests, so that it obtains varius
# KDEs for the positions of enzymes, and uses those KDEs to smooth the curve and reduce the variability
# on the correlation.

# Inputs
# --------------------------------------------------------------
file_out = 'test_calibration_nsets'  # p_calibraiton =processed calibration

# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 24#60#120#24#12#120#60  #12#24  #24#8 #96 # 120
n_sets = 20

# Simulation conditions
# --------------------------------------------------------------
dt = 0.25
initial_time = 0
final_time = 1000#1000#2000#2000#500#2000#500 # o 2000
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Reference - It is the reference density of topos when there is no gene that we will use to calculate the
#             fold enrichment.
# --------------------------------------------------------------
reference_path = '../noRNAP/'  # Path with the reference kdes

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
colors_dict = {'topoI': 'red', 'gyrase': 'cyan'}
colors_dict = {'RNAP': 'purple', 'topoI': 'red', 'gyrase': 'cyan'}
kwargs = {'linewidth': 2, 'ls': '-'}
ylabel = 'Fold Enrichment'
ylabel2 = 'RNAP Density'
#ylabels = ['Density', 'Fold-Enrichment', 'Fold-Enrichment']

# TARGETS FOR OPTIMIZATION - We need them for hte parallelization
# -----------------------------------
target_FE = 1.68  # Target fold-enrichment
target_CO = 1.0  # Target correlation between topo I and RNAP densities.
x_spacing = 10.0  # The spacing I want at the moment of doing interpolation.

# nbins is the number of bins to use when calculating the kde
target_dict = {'target_FE': target_FE, 'target_CO': target_CO, 'target_gene': 'reporter',
               'enzymes_names': ['RNAP', 'topoI', 'gyrase']}


# ----------------------------------------------------------------------------------------------------------------------
# Parallelization functions
# ----------------------------------------------------------------------------------------------------------------------

# This time we call it parallelization_function to distinguish from the calibration/optimizaiton function.
def parallelization_function():
    # We need to prepare the inputs.
    # At the moment, we only have one system.

    # Global dictionaries
    # ------------------------------------------
    global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                   'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                   'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                   'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'n_sets': n_sets,
                   # 'initial_sigma': initial_sigma,
                   'DNA_concentration': 0.0}

    # Create lists of conditions for each system
    # ------------------------------------------

    # Global dictionaries
    global_dict_list = [global_dict]

    # List of lists of variations
    variations_list = [[]]

    # Arrays with position densities to calculate fold change
    list_reference = [reference_dict]

    # Finally, run objective function.
    # ------------------------------------------
    my_objective, output_dict = tct.single_case_RNAPTracking_calibration_nsets(global_dict_list, variations_list,
                                                                               list_reference, n_simulations,
                                                                               target_dict)

    return my_objective, output_dict


# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------

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

reference_dict = {}
for name in target_dict['enzymes_names']:

    # Load reference file
    if name != 'RNAP':
        kde_ref = np.loadtxt(reference_path + 'reference_' + name + '.txt')

        # And get interpolated data - I do it here, so we don't have to do it again in the parallelization and save some
        # time. It is saved to the reference_dict by the way
        reference_dict[name] = tct.get_interpolated_kde(kde_ref[:, 0], kde_ref[:, 1], x_system)

# Run parallelization
# -----------------------------------------------------
objective, output = parallelization_function()

# Let's plot the densities kdes and histograms
# ---------------------------------------------------------
xs = np.arange(initial_time, final_time + 2 * dt, dt)  # For plotting superhelical

kde_overall = {}
fig, axs = plt.subplots(4, figsize=(width, 4 * height), tight_layout=True)
# 4 plots: RNAP, topoI, gyrase and superhelical density

for p, name in enumerate(names):

    ax = axs[p]

    # Unpack histograms
    counts = output['results']['hists_dict'][name]['mean']
    counts_std = output['results']['hists_dict'][name]['std']
    bin_edges = output['results']['hists_dict'][name]['bin_edges']

    # Calculate the width of each bin
    bin_width = bin_edges[1] - bin_edges[0]

    # Plot bars for each bin
    ax.bar(bin_edges[:-1], counts, width=bin_width, color=colors_dict[name], edgecolor='black', alpha=0.5)

    # Different process for RNAP or topos
    if name == 'RNAP':
        x = x_gene
        y = output['results']['kde_gene'][name]['mean']
        ys = output['results']['kde_gene'][name]['std']
    else:
        x = x_system
        y = output['results']['kde_system'][name]['mean']
        ys = output['results']['kde_system'][name]['std']

    ax.plot(x, y, colors_dict[name], lw = lw, label = name)
    ax.fill_between(x, y-ys, y+ys, color=colors_dict[name], alpha=0.5)
    #if name == 'topoI':
    #    ax2.plot(x_system, output['results']['kde_system'][name], colors_dict[name], lw=lw, label=name)

    # Labels
    # ------------------------
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    #ax.set_title(name, fontsize=title_size)
    #ax.set_ylim(0,np.max(kde_y))
    ax.set_xlim(0, my_circuit.size)
    ax.grid(True)
    ax.legend(loc='best', fontsize=font_size)
    #ax2.legend(loc='best', fontsize=font_size)


print('Averaged FE:', output['results']['FE_val']['topoI']['mean'], output['results']['FE_val']['topoI']['std'])
print('Correlation:', output['results']['correlation']['mean'], output['results']['correlation']['std'])
print('Overall correlation', output['overall_correlation'])
print('objective', objective)

# And plot superhelicals
# ------------------------------------------------
ax = axs[3]
y = output['results']['superhelical_dict']['mean'][:]
ys = output['results']['superhelical_dict']['std'][:]
y1 = y - ys
y2 = y + ys

ax.plot(xs, y, 'black', lw=lw)
ax.fill_between(xs, y1, y2, color='black', alpha=0.5)
ax.set_ylabel(r'Global $\sigma$', fontsize=font_size)
ax.set_xlabel(r'Time (seconds)', fontsize=font_size)
ax.grid(True)

plt.savefig(file_out + '.png')

plt.show()
