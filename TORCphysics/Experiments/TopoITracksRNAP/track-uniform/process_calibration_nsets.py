import numpy as np
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import matplotlib.pyplot as plt


# Description
# --------------------------------------------------------------
# We want to plot results to check that the graph we want to produce at the end makes sense, using the
# calibration using nsets to smooth the KDEs

# Inputs
# --------------------------------------------------------------
file_out = 'pcalibration_nsets'  # p_calibraiton =processed calibration

# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 12#120#60  #12#24  #24#8 #96 # 120
n_sets = 20#0

# Simulation conditions
# --------------------------------------------------------------
dt = 0.25
initial_time = 0
final_time = 500 # o 2000
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
                   'frames': frames, 'dt': dt, 'n_simulations': n_simulations,  'n_sets': n_sets, # 'initial_sigma': initial_sigma,
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

    #if name == 'topoI':
    #    ax2.plot(x_system, output['results']['kde_system'][name], colors_dict[name], lw=lw, label=name)

# Labels
# ------------------------
ax.set_ylabel(ylabel, fontsize=font_size)
ax2.set_ylabel(ylabel2, fontsize=font_size, color=colors_dict['RNAP'])
ax.set_xlabel(r'Position (bp)', fontsize=font_size)
#ax.set_title(name, fontsize=title_size)
#ax.set_ylim(0,np.max(kde_y))
ax.set_xlim(0, my_circuit.size)
ax.grid(True)
ax.legend(loc='best', fontsize=font_size)
#ax2.legend(loc='best', fontsize=font_size)

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
plt.savefig(file_out + '.pdf')
plt.savefig(file_out + '.png')
print('FE',FE)
print('CO',CO)
print('objective', objective)
plt.show()
