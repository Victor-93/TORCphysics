import numpy as np
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import matplotlib.pyplot as plt
import pickle
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns
import sys

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
percentage_threshold = .005

# Inputs
# --------------------------------------------------------------

pickle_files = [
    # Second avg test - with Stages-Stall fixed
    #'track-StagesStall/avg02/calibration_avg-RNAPTracking_nsets_p2_small_dt'+str(dt)+'.pkl'
#    'track-StagesStall/avg02/avgx2-reproduce-calibration_RNAPTracking_nsets_p2_small_dt'+str(dt)+'.pkl'
    #'track-StagesStall/small_distance/02-calibration_avg-RNAPTracking_nsets_p2_small_dt' + str(dt) + '.pkl'
#    'track-StagesStall/test-borrar-trials'+'.pkl'
    [f'track-StagesStall/trials_bigdist_v2/0{i}-02-calibration_avg-RNAPTracking_nsets_p2_small_dt1.0-trials.pkl' for i in range(1, 3)],
    [f'track-StagesStall/trials_bigdist/0{i}-calibration_avg-RNAPTracking_nsets_p2_small_dt1.0-trials.pkl' for i in range(1, 5)],
    [f'track-StagesStall/trials_smalldist/0{i}-calibration_avg-RNAPTracking_nsets_p2_small_dt1.0_smalldist-trials.pkl' for i in range(1, 5)]
]

system_labels = ['bid_dist_v2', 'bid_dist', 'small_dist']

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
fig, axs = plt.subplots(6,ncases, figsize=(ncases*width, 6*height), tight_layout=True)
outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']



x1 = 3500
x0 = 2500#-30
h = [2.3, 2.3, 1.7]
dx = x1 - x0
gene_colour = 'gray'
gene_lw=3

ylims = [[.7,2.5], [.7,2.5], [.67,1.8]]
#ylims = [[.75,1.8]]
#ylims = [[.0,1]]
mylabels = ['RNAP', 'Topoisomerase I', 'Gyrase']

for k, pfiles in enumerate(pickle_files):


    # Load data and process losses according threshold
    # ---------------------------------------------------------
    trials_data = []
    for pfile in pfiles:
        with open(pfile, 'rb') as file:
            trials_data.extend( pickle.load(file))

    results_data = []
    for pfile in pfiles:
        with open(pfile, 'rb') as file:
            results_data.extend( pickle.load(file).results)

    # Loss distributions and filtered percentage
    # ---------------------------------------------------------
    # Let's sort the losses
    results = results_data

    loss_df = pd.DataFrame({'loss': [t['loss'] for t in results]})

    loss_df = loss_df.sort_values(by='loss', ascending=False)  # , inplace=True)

    n = len(loss_df['loss'])
    nconsidered = int(n * percentage_threshold)
    err_threshold = loss_df['loss'].iloc[-nconsidered]
    print('Number of tests', n)
    print('Considered', nconsidered)
    print('For ', percentage_threshold * 100, '%')
    # Filter according error
    filtered_df = loss_df[loss_df['loss'] <= err_threshold]

    # Set with minimum loss
    dat = min(results, key=lambda x: x['loss'])

    # Let's filter the range of values that give the best result
    filtered_trials = [trial for trial in trials_data if trial['result']['loss'] <= err_threshold] # This contains all info
    # Now, the dict with the actual parametrisation from the random search
    filtered_oparams_dict = [
        {key: val[0] for key, val in trial['misc']['vals'].items()}
        for trial in filtered_trials
    ]

    # filtered_oparams_dict = [trial['misc']['vals'] for trial in filtered_trials]
    filtered_oparams_df = pd.DataFrame(filtered_oparams_dict)  # And make it a dataframe

    # And the best trial
    best_trial = min(filtered_trials, key=lambda x: x['result']['loss'])

    # Loss distributions and filtered percentage
    # ---------------------------------------------------------
    ax = axs[0,k]
    ax.set_title('Loss distribution ' + system_labels[k])

    sns.violinplot(data=loss_df, ax=ax, inner="quart")  # , cut=0, color=colors[i])
 #   ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel('Loss')
    ax.grid(True)

    ax.set_ylim([0,4])
    #ax.set_yscale('log')

    # Filtered loss distribution
    # ----------------------------------------------------------------------------------------------------------------------
    ax = axs[1,k]
    ax.set_title('Filtered loss for ' + system_labels[k])

    sns.violinplot(data=filtered_df, ax=ax, inner="quart")  # , cut=0, color=colors[i])
#    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel('Loss')
    ax.grid(True)

    # Histogram
    # ----------------------------------------------------------------------------------------------------------------------
    ax = axs[2,k]
    loss = loss_df['loss'].to_numpy()
    floss = filtered_df['loss'].to_numpy()

    # Create a histogram
    minv = min(loss)
    maxv = np.mean(loss) + 1 * np.std(loss)
    maxv = 0.075#.5  # max(loss)*.2
    bins = np.linspace(minv, maxv, 100)  # Define bins
    hist, bin_edges = np.histogram(loss, bins=bins)

    # Plot the full histogram
    ax.hist(loss, bins=bins, color='gray', alpha=0.6, label='Loss')

    # Highlight bins corresponding to floss
    for value in floss:
        # Find the bin index for the current value
        bin_index = np.digitize(value, bin_edges) - 1
        # Plot the specific bin
        ax.bar(
            bin_edges[bin_index],  # Bin start
            hist[bin_index],  # Bin height
            width=bin_edges[1] - bin_edges[0],  # Bin width
            color='red',  # Highlight color
            alpha=0.8,
            #        edgecolor='black',
            label='Highlighted' if bin_index == np.digitize(floss[0], bin_edges) - 1 else ""
        )
    # ax.set_xlabel('test')
    # ax.set_ylabel('loss')
    #ax.set_xscale('log')
    ax.grid(True)

    # ----------------------------------------------------------------------------------------------------------------------
    # Best result
    # ----------------------------------------------------------------------------------------------------------------------
    ax = axs[3,k]

    # Draw gene
    # ---------------------------------------------------------
    arrow = mpatches.FancyArrowPatch((x0, h[k]), (x1, h[k]), arrowstyle='simple',
                                     facecolor=gene_colour, zorder=5, edgecolor='black', lw=gene_lw,
                                     mutation_scale=40, shrinkA=0, shrinkB=0)
    ax.add_patch(arrow)

    output = dat #min(data, key=lambda x: x['loss']) # Select the one with the minimum loss

    for p, name in enumerate(names):

        ax = axs[3, k]

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
    yl = ylims[k]
    ax.set_ylim(yl[0],yl[1])
    ax.set_xlim(0, my_circuit.size)
    #ax.set_title(title[n], fontsize=title_size)
    ax2.tick_params(axis='y', labelcolor=colors_dict['RNAP'])
    ax2.set_ylim(0,1)

    # fig.suptitle(title[n], fontsize=title_size)
    if k == 0:
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
#    print('For ', title[n])
    print('FE',FE)
    print('CO',CO)
    print('objective', output['loss'])
    print('RNAP_correlation', output['RNAP_correlation'])

    ax.plot([x0,x0], [-10,10], lw=2, color='gray', ls='--' )
    ax.plot([x1, x1], [-10, 10], lw=2, color='gray', ls='--')

    # Add label outside the plot
    ax.text(-0.1, 1.1, outside_label[0], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

    # ----------------------------------------------------------------------------------------------------------------------
    # Values - best set
    # ----------------------------------------------------------------------------------------------------------------------
    ax = axs[4,k]
    ax.set_title('Best set vals')

    case_dict = filtered_oparams_df #pd.read_csv(param_file)
    case_dict['RNAP_dist'] = case_dict['RNAP_dist'].apply(lambda x: x/1000)
    case_dict['fold_change'] = case_dict['fold_change'].apply(lambda x: x / 100)
    mean_vals = case_dict.mean(axis=0)
    std_vals = case_dict.std(axis=0)

    # Overalls
    ov_df = pd.concat([mean_vals, std_vals], axis=1)

    print('mean vals:', ov_df)
#    case_dict['RNAP_dist'] =case_dict['RNAP']/100.0
 #   nokeys = ['RNAP_dist', 'fold_change', 'gamma']

    # Remove the keys from all dictionaries in the list
#    for key in nokeys:
#        if key in case_dict:
#            del case_dict[key]

    keys = list(case_dict.keys())

    # Number of cases and keys
    n_cases = 1
    n_keys = len(keys)

    # Set up the positions for each bar group (x-axis)
    x = np.arange(n_keys)  # Position of each group on the x-axis
    width = 0.8 / n_cases  # Dynamically calculate the width of each bar

#    keys = ['spacer_kon', 'k_off', 'k_open', 'k_closed', 'k_ini']

    # Plot each case in the list
    values_m = np.array([mean_vals[key] for key in keys])  # Get values for the current case
    values_s = np.array([std_vals[key] for key in keys])  # Get values for the current case
    ax.bar(x +  width / 2, values_m, width, yerr=values_s, color='purple')

    # Add labels, title, and custom ticks
    ax.set_xlabel('Rate Type', fontsize=xlabel_size)
    ax.set_ylabel(r'Rate ($s^{-1}$)', fontsize=xlabel_size)
    ax.set_title('Transition Rates - Best set', fontsize=title_size)
    ax.set_xticks(x)
    ax.set_xticklabels(keys)
    ax.grid(True)

    # Add label outside the plot
    ax.text(-0.1, 1.1, outside_label[1], transform=ax.transAxes,
            fontsize=font_size * 1.5, fontweight='bold', va='center', ha='center')

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="center")

    #axs[2,1].text(-0.1, 1.1, outside_label[5], transform=axs[2,1].transAxes,
    #        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

    # ----------------------------------------------------------------------------------------------------------------------
    # Values - best trial
    # ----------------------------------------------------------------------------------------------------------------------
    ax = axs[5,k]
    ax.set_title('Best set vals')

    #case_dict = [
    #    {key: val[0] for key, val in [best_trial].items()}
    #
    #for trial in [best_trial]
   # ]
    case_dict = best_trial['misc']['vals']

    for val in case_dict:
        case_dict[val] = case_dict[val][0]
        if 'RNAP_dist' in val:
            case_dict[val] = case_dict[val]/1000
        if 'fold_change' in val:
            case_dict[val] = case_dict[val]/100

    keys = list(case_dict.keys())

    # Number of cases and keys
    n_cases = 1
    n_keys = len(keys)

    # Set up the positions for each bar group (x-axis)
    x = np.arange(n_keys)  # Position of each group on the x-axis
    width = 0.8 / n_cases  # Dynamically calculate the width of each bar

    # Plot each case in the list
    values = np.array([case_dict[key] for key in keys])  # Get values for the current case
    ax.bar(x +  width / 2, values, width, color='red')

    # Add labels, title, and custom ticks
    ax.set_xlabel('Rate Type', fontsize=xlabel_size)
    ax.set_ylabel(r'Rate ($s^{-1}$)', fontsize=xlabel_size)
    ax.set_title('Transition Rates - Best case', fontsize=title_size)
    ax.set_xticks(x)
    ax.set_xticklabels(keys)
    ax.grid(True)

    # Add label outside the plot
    ax.text(-0.1, 1.1, outside_label[1], transform=ax.transAxes,
            fontsize=font_size * 1.5, fontweight='bold', va='center', ha='center')

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="center")

#plt.savefig(output_prefix+'-FE.png')
#plt.savefig(output_prefix+'-FE.pdf')


plt.show()
