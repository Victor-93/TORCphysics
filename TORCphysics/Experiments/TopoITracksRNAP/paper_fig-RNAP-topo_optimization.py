import numpy as np
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import matplotlib.pyplot as plt
import pickle
import matplotlib.patches as mpatches

# Description
# --------------------------------------------------------------
# Following the calibration process by executing calibrate_tracking_nsets_p2.py for both the stalling and uniform
# models, here I plot their results

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
    'track-uniform/calibration_RNAPTracking_nsets_p2_small_dt'+str(dt)+'.pkl',
    'track-StagesStall/reproduce-calibration_RNAPTracking_nsets_p2_small_dt'+str(dt)+'-02.pkl'
]
output_prefix = 'RNAPStages-topoIRNAPtracking'
title = ['Uniform Model', 'Stages-Stall Model']

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
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

alpha=0.3

names = ['RNAP', 'topoI', 'gyrase']
colors_dict = {'RNAP': 'purple', 'topoI': 'red', 'gyrase': 'blue'}
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

x1 = 3500
x0 = 2500#-30
h = [.725,.725]
dx = x1 - x0
gene_colour = 'gray'
gene_lw=2

ylims = [[.67,1.5], [.67,1.8]]
mylabels = ['RNAP', 'Topoisomerase I', 'Gyrase']

for n, pickle_file in enumerate(pickle_files):

    # Draw gene
    # ---------------------------------------------------------
    arrow = mpatches.FancyArrowPatch((x0, h[n]), (x1, h[n]), arrowstyle='simple',
                                     facecolor=gene_colour, zorder=5, edgecolor='black', lw=gene_lw,
                                     mutation_scale=40, shrinkA=0, shrinkB=0)
    axs[n].add_patch(arrow)

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
            ymax = np.max(y)
            ys = output['results']['kde_gene'][name]['std']
            y=y/ymax
            ys=ys/ymax
            # ax2 = ax.twinx()
            # ax2.plot(x_gene, output['results']['kde_gene'][name]['mean'], colors_dict[name], lw=lw, label=name)
        else:
            x = x_system
            y = output['results']['FE_curve'][name]['mean']
            ys = output['results']['FE_curve'][name]['std']
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
    FE = output['FE'][0]
    CO = output['overall_correlation']
    OB = output['objective']
    RNAP_CO = output['results']['RNAP_correlation']

    # Define the text to display
    # textstr = f'FE={FE:.2f}, RCO={RNAP_CO:.2f}, CO={CO:.2f}, OB={OB:.2e}'
    textstr = f'FE={FE:.2f}, RCO={RNAP_CO:.2f}, CO={CO:.2f}' #, OB={OB:.2e}'

    # Add the text box to the plot
    props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
    ax.text(0.6, .93, textstr, transform=ax.transAxes, fontsize=font_size*.95,
            verticalalignment='top', bbox=props)
    print('For ', title[n])
    print('FE',FE)
    print('CO',CO)
    print('objective', output['objective'])
    print('RNAP_correlation', output['results']['RNAP_correlation'])

    ax.plot([x0,x0], [-10,10], lw=2, color='gray', ls='--' )
    ax.plot([x1, x1], [-10, 10], lw=2, color='gray', ls='--')
plt.savefig(output_prefix+'-FE.png')
plt.savefig(output_prefix+'-FE.pdf')
plt.show()
