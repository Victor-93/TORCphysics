import numpy as np
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import matplotlib.pyplot as plt


# Description
# --------------------------------------------------------------
# This is a script that is very similar to process_calibration.py, but here I just want to make sure
# that my calibration is working and the parameters I am using to optimize are sensible and deterministic enough
# so our calibration process is trustworthy. This because at the moment, I get very variable correlation values.

# TODO: Aqui testeemos cosas.
#  1.- Calcula avg FE en la gene region
#  2.- Outputea del calibrationtools additional_results = {}, que contenga supercoiling, y histogram.
#  3.- Checa por que da tan feo el kde del RNAP.
#  4.- Outputea los CO del hyperopt

# TODO: De ahorita para asegurarnos que sirve nuestro proceso.
#  * Plotea histograma y kde para ver qu estamos calculando lo que queremos correctamente.
#  * Varia el tiempo y el numero de simulaciones para que veas mas o menos cuanto necesitas para obtener distribuciones
#    similares.
# Inputs
# --------------------------------------------------------------
file_out = 'test_calibration'  # p_calibraiton =processed calibration

# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 36#60#120#24#12#120#60  #12#24  #24#8 #96 # 120
n_tests = 10

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
                   'frames': frames, 'dt': dt, 'n_simulations': n_simulations,  # 'initial_sigma': initial_sigma,
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
    my_objective, output_dict = tct.single_case_RNAPTracking_calibration(global_dict_list, variations_list,
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
RNAP_env =  [environment for environment in my_circuit.environmental_list if environment.name == 'RNAP'][0]

# Define x-axes
x_system = tct.get_interpolated_x(1, my_circuit.size)
#x_gene = tct.get_interpolated_x(target_gene.start, target_gene.end)
x_gene = tct.get_interpolated_x(target_gene.start-RNAP_env.size, target_gene.end)

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
objective_list = []
output_list = []
for n in range(n_tests):
    objective, output = parallelization_function()
    objective_list.append(objective)
    output_list.append(output)

# Let's plot the densities kdes and histograms
# ---------------------------------------------------------
#fig, axs = plt.subplots(3,len(n_tests), figsize=(width, 1*height), tight_layout=True)
xs = np.arange(initial_time, final_time + 2*dt, dt) # For plotting superhelical

kde_overall = {}
fig, axs = plt.subplots(4, n_tests+1, figsize=(width*(n_tests+1), 4*height), tight_layout=True)
for n in range(n_tests):
    output = output_list[n]
    objective = objective_list[n]

    for p, name in enumerate(names):

        if n == 0:
            kde_overall[name] = np.zeros_like(output['results']['kde_gene'][name])

        kde_overall[name] = kde_overall[name] + output['results']['kde_gene'][name]

        ax = axs[p,n]

        # Unpack histograms
        counts = output['results']['hists_dict'][name]['counts']
        bin_edges = output['results']['hists_dict'][name]['bin_edges']

        # Calculate the width of each bin
        bin_width = bin_edges[1] - bin_edges[0]

        # Plot bars for each bin
        ax.bar(bin_edges[:-1], counts, width=bin_width, color=colors_dict[name], edgecolor='black', alpha=0.5)

        # Different process for RNAP or topos
        if name == 'RNAP':
            ax.plot(x_gene, output['results']['kde_gene'][name], colors_dict[name], lw=lw, label=name)
        else:
            ax.plot(x_system, output['results']['kde_system'][name], colors_dict[name], lw=lw, label=name)

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

        # Add FE and correlation labels
        # -------------------------------------
        FE = output['FE']
        CO = output['correlation']

        # Define the text to display
        textstr = f'FE={FE:.2f}, CO={CO:.2f}'

        # Add the text box to the plot
        props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
        ax.text(0.05, 0.1, textstr, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
    #plt.savefig(file_out + '.pdf')
    #plt.savefig(file_out + '.png')
    print('FE',FE)
    print('CO',CO)
    print('objective', objective)

    # And plot superhelicals
    # ------------------------------------------------
    ax = axs[3,n]
    y = output['results']['superhelical_dict']['mean'][:]
    ys = output['results']['superhelical_dict']['std'][:]
    y1 = y-ys
    y2 = y+ys

    print(len(xs), len(y1), len(y2))

    ax.plot(xs,y, 'black', lw=lw)
    ax.fill_between(xs, y1, y2, color='black', alpha=0.5)
    ax.set_ylabel(r'Global $\sigma$', fontsize=font_size)
    ax.set_xlabel(r'Time (seconds)', fontsize=font_size)
    ax.grid(True)


# And plot kde overalls ( we need to smooth the curves)
# ------------------------------------------------
for p, name in enumerate(names):

    print(p,n_tests+1)

    ax = axs[p, n_tests]

    kde_overall[name] = kde_overall[name]/n_tests
    buf_size = 10

    # Different process for RNAP or topos
    ax.plot(x_gene, kde_overall[name], colors_dict[name], lw=lw, label=name)

    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    ax.set_xlim(0, my_circuit.size)
    ax.grid(True)

ax = axs[0,n_tests]
correlation_matrix = np.corrcoef(kde_overall['topoI'][buf_size:-buf_size], kde_overall['RNAP'][buf_size:-buf_size])
correlation = correlation_matrix[0, 1]

# Define the text to display
textstr = f'OvCO={correlation:.2f}'

# Add the text box to the plot
props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
ax.text(0.05, 0.1, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

print('overall correlation:', correlation)
plt.savefig(file_out + '.pdf')


plt.show()
