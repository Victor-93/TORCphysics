import numpy as np
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import matplotlib.pyplot as plt
import pickle
import matplotlib.patches as mpatches
import pandas as pd

# Description
# --------------------------------------------------------------
# Following the calibration process by executing calibrate_tracking_nsets_p2.py for both the stalling and uniform
# models, here I plot their results

# TODO: You have to find a way to figure out the averaged values of rates.

# Simulation conditions
# --------------------------------------------------------------
dt =  1.0
initial_time = 0
final_time = 200 # 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Inputs
# --------------------------------------------------------------

pickle_files = [
    # Second avg test - with Stages-Stall fixed
    #'track-StagesStall/avg02/calibration_avg-RNAPTracking_nsets_p2_small_dt'+str(dt)+'.pkl'
#    'track-StagesStall/avg02/avgx2-reproduce-calibration_RNAPTracking_nsets_p2_small_dt'+str(dt)+'.pkl'
    #'track-StagesStall/small_distance/02-calibration_avg-RNAPTracking_nsets_p2_small_dt' + str(dt) + '.pkl'
    'track-StagesStall/test-borrar-trials'+'.pkl'
]
param_file = 'avgx2_small-dist_table_dt1.0.csv'

output_prefix = 'RNAPStages-topoIRNAPtracking_v2_wrates'
title = ['Topoisomerase I follows RNAP experiment']

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
height = 3.5
lw = 2
font_size = 12
xlabel_size = 14
title_size = 16

alpha=0.3

names = ['RNAP', 'topoI', 'gyrase']
colors_dict = {'RNAP': 'purple', 'topoI': 'red', 'gyrase': 'blue'}
kwargs = {'linewidth': 2, 'ls': '-'}
ylabel = 'Topo Fold Enrichment'
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
fig, axs = plt.subplots(ncases+1, figsize=(width, 2*height), tight_layout=True)
outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

x1 = 3500
x0 = 2500#-30
h = [2.3]
dx = x1 - x0
gene_colour = 'gray'
gene_lw=3

#ylims = [[.67,1.5], [.67,1.8]]
ylims = [[.75,1.8]]
#ylims = [[.0,1]]
mylabels = ['RNAP', 'Topoisomerase I', 'Gyrase']

for n, pickle_file in enumerate(pickle_files):

    # Draw gene
    # ---------------------------------------------------------
    arrow = mpatches.FancyArrowPatch((x0, h[n]), (x1, h[n]), arrowstyle='simple',
                                     facecolor=gene_colour, zorder=5, edgecolor='black', lw=gene_lw,
                                     mutation_scale=40, shrinkA=0, shrinkB=0)
    axs[0].add_patch(arrow)

    # Load the dictionary from the file
    # ---------------------------------------------------------
    with open(pickle_file, 'rb') as file:
        #output = pickle.load(file)
        data = pickle.load(file).results

    # Filter data to include only entries with 'status' == 'ok'
    data = [item for item in data if item['status'] == 'ok']

    output = min(data, key=lambda x: x['loss']) # Select the one with the minimum loss

    for p, name in enumerate(names):

        ax = axs[0]

        # Different process for RNAP or topos
        if name == 'RNAP':
            ax2 = ax.twinx()  # Different scale for RNAP
            ax = ax2
            x = x_gene
            y = output['kde_gene'][name]['mean']
            ymax = np.max(y)
            ys = output['kde_gene'][name]['std']
            y=y/ymax
            ys=ys/ymax
            # ax2 = ax.twinx()
            # ax2.plot(x_gene, output['results']['kde_gene'][name]['mean'], colors_dict[name], lw=lw, label=name)
        else:
            x = x_system
            y = output['FE_curve'][name]['mean']
            ys = output['FE_curve'][name]['std']
            # ax.plot(x_system, output['results']['FE_curve'][name]['mean'], color=colors_dict[name], lw=lw, label=name)
        ax.plot(x, y, color=colors_dict[name], lw=lw, label=mylabels[p])
        ax.fill_between(x, y-ys, y+ys, color=colors_dict[name], alpha=alpha)
        print(name, np.mean(y), p)

    # Labels
    # ------------------------
    ax.set_ylabel(ylabel, fontsize=xlabel_size)
    ax2.set_ylabel(ylabel2, fontsize=xlabel_size, color=colors_dict['RNAP'])
    ax.set_xlabel(r'Position (bp)', fontsize=xlabel_size)
    ax.set_xlim(0, my_circuit.size)
    ax.grid(True)
    yl = ylims[n]
    ax.set_ylim(yl[0],yl[1])
    ax.set_xlim(0, my_circuit.size)
    ax.set_title(title[n], fontsize=title_size)
    ax2.tick_params(axis='y', labelcolor=colors_dict['RNAP'])
    ax2.set_ylim(0,1)

    # fig.suptitle(title[n], fontsize=title_size)
    if n == 0:
        # Combine the legends from both axes
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()

        # Create one unified legend
        ax.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=font_size)
        #ax.legend(loc='best', fontsize=font_size)
        #ax2.legend(loc='best', fontsize=font_size)


    # Add FE and correlation labels
    # -------------------------------------
    FE = output['FE'] # [0]
    CO = output['overall_correlation']
    # OB = output['objective']
    OB = output['loss']
    RNAP_CO = output['RNAP_correlation']

    # Define the text to display
    # textstr = f'FE={FE:.2f}, RCO={RNAP_CO:.2f}, CO={CO:.2f}, OB={OB:.2e}'
    # textstr = f'FE={FE:.2f}, RCO={RNAP_CO:.2f}, CO={CO:.2f}' #, OB={OB:.2e}'
    # textstr = f'$\\mathcal{{F}}={FE:.2f}$, $\\lambda={RNAP_CO:.2f}$, $\\rho={CO:.2f}$'
    textstr = f'$\\mathcal{{F}}={FE:.2f}$, $\\rho={CO:.2f}$, $\\lambda={RNAP_CO:.2f}$ \n $f={OB:.2e}$'

    # Add the text box to the plot
    props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
    ax.text(0.63, .93, textstr, transform=ax.transAxes, fontsize=font_size*.95,
            verticalalignment='top', bbox=props)
    print('For ', title[n])
    print('FE',FE)
    print('CO',CO)
    print('objective', output['loss'])
    print('RNAP_correlation', output['RNAP_correlation'])

    ax.plot([x0,x0], [-10,10], lw=2, color='gray', ls='--' )
    ax.plot([x1, x1], [-10, 10], lw=2, color='gray', ls='--')

    # Add label outside the plot
    ax.text(-0.1, 1.1, outside_label[0], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')


# Rates
ax = axs[1]
case_dict = pd.read_csv(param_file)
nokeys = ['RNAP_dist', 'fold_change', 'gamma']

# Remove the keys from all dictionaries in the list
#for case_dict in responses:
for key in nokeys:
    if key in case_dict:
        del case_dict[key]

keys = list(case_dict.keys())

# Number of cases and keys
n_cases = 1
n_keys = len(keys)

# Set up the positions for each bar group (x-axis)
x = np.arange(n_keys)  # Position of each group on the x-axis
width = 0.8 / n_cases  # Dynamically calculate the width of each bar

keys = ['k_on', 'k_off', 'k_open', 'k_closed', 'k_ini']

# Plot each case in the list
values = np.array([case_dict[key] for key in keys])  # Get values for the current case
ax.bar(x +  width / 2, values[:,0], width, yerr=values[:,1], color='purple')

# Add labels, title, and custom ticks
ax.set_xlabel('Rate Type', fontsize=xlabel_size)
ax.set_ylabel(r'Rate ($s^{-1}$)', fontsize=xlabel_size)
ax.set_title('Transition Rates', fontsize=title_size)
ax.set_xticks(x)
ax.set_xticklabels(keys)
ax.grid(True)

# Add label outside the plot
ax.text(-0.1, 1.1, outside_label[1], transform=ax.transAxes,
        fontsize=font_size * 1.5, fontweight='bold', va='center', ha='center')

#axs[2,1].text(-0.1, 1.1, outside_label[5], transform=axs[2,1].transAxes,
#        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')


#plt.savefig(output_prefix+'-FE.png')
#plt.savefig(output_prefix+'-FE.pdf')


plt.show()
