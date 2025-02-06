import numpy as np
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import matplotlib.pyplot as plt
import pickle
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns
import sys
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.cluster import KMeans

# Description
# --------------------------------------------------------------
# I will analyse trials results, but this time I will do a kmeans to see if I can get better statistics, as it seems
# that using the loss as the only parameter does not give good enough results. This is because maybe there are
# multiple solutions.

# Simulation conditions
# --------------------------------------------------------------
dt =  1.0
initial_time = 0
final_time = 200 # 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
percentage_threshold = 0.05 #This time is bigger because we want to get rid of half of worst results
n_clusters = [3,2]

# Inputs
# --------------------------------------------------------------

pickle_files = [
    # Second avg test - with Stages-Stall fixed
    #'track-StagesStall/avg02/calibration_avg-RNAPTracking_nsets_p2_small_dt'+str(dt)+'.pkl'
#    'track-StagesStall/avg02/avgx2-reproduce-calibration_RNAPTracking_nsets_p2_small_dt'+str(dt)+'.pkl'
    #'track-StagesStall/small_distance/02-calibration_avg-RNAPTracking_nsets_p2_small_dt' + str(dt) + '.pkl'
#    'track-StagesStall/test-borrar-trials'+'.pkl'
    [f'track-StagesStall/trials_bigdist/0{i}-calibration_avg-RNAPTracking_nsets_p2_small_dt1.0-trials.pkl' for i in range(1, 5)],
    [f'track-StagesStall/trials_smalldist/0{i}-calibration_avg-RNAPTracking_nsets_p2_small_dt1.0_smalldist-trials.pkl' for i in range(1, 5)]
]

system_labels = ['bid_dist', 'small_dist']

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
fig, axs = plt.subplots(2,ncases, figsize=(ncases*width, 2*height), tight_layout=True)
outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']



x1 = 3500
x0 = 2500#-30
h = [2.3, 1.7]
dx = x1 - x0
gene_colour = 'gray'
gene_lw=3

ylims = [[.7,2.5], [.67,1.8]]
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
    err_threshold = loss_df['loss'].iloc[-nconsidered] # This is the error threshold

    # Let's filter the range of values that give the best result
    filtered_trials = [trial for trial in trials_data if trial['result']['loss'] <= err_threshold] # This contains all info
    filtered_trials_loss = [trial['result']['loss'] for trial in filtered_trials]

    # Now, the dict with the actual parametrisation from the random search
    filtered_oparams_dict = [
        {key: val[0] for key, val in trial['misc']['vals'].items()}
        for trial in filtered_trials
    ]

    # filtered_oparams_dict = [trial['misc']['vals'] for trial in filtered_trials]
    filtered_oparams_df = pd.DataFrame(filtered_oparams_dict)  # And make it a dataframe

    # Normalise the data
    # --------------------------------------------------------------------------------------
    scaler = MinMaxScaler()  # Use MinMaxScaler() if you prefer scaling to [0, 1]
    normalized_data = scaler.fit_transform(filtered_oparams_df)

    # Convert the normalised data back to a DataFrame for easy inspection (optional)
    normalized_df = pd.DataFrame(normalized_data, columns=filtered_oparams_df.columns)

    # Clustering
    # --------------------------------------------------------------------------------------
    # Set the number of clusters
#    n_clusters = 3  # Change this to the desired number of clusters

    # Apply K-Means
    kmeans = KMeans(n_clusters=n_clusters[k], random_state=42)
    kmeans.fit(normalized_data)

    # Get cluster labels
    filtered_oparams_df['cluster'] = kmeans.labels_
    filtered_oparams_df['loss'] = filtered_trials_loss

    # Ensure results are fully calculated
    assert 'cluster' in filtered_oparams_df.columns, "Clustering not complete!"

    # Inspect the clusters
#    print(filtered_oparams_df.head())

    # Calculate mean and std for each cluster
    cluster_summary = filtered_oparams_df.groupby('cluster').mean()
    print(cluster_summary)

    #print('Number of tests', n)
    #print('Considered', nconsidered)
    #print('For ', percentage_threshold * 100, '%')
    # Filter according error
    #filtered_df = loss_df[loss_df['loss'] <= err_threshold]

    # Set with minimum loss
    #dat = min(results, key=lambda x: x['loss'])

    # And the best trial
    #best_trial = min(filtered_trials, key=lambda x: x['result']['loss'])

    # Loss distributions and filtered percentage
    # ---------------------------------------------------------
    ax = axs[0,k]
    ax.set_title('Loss distribution ' + system_labels[k])

    # Create violin plot
    sns.violinplot(data=filtered_oparams_df, hue='cluster', y='loss', inner="quart", palette="viridis",ax=ax)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel('Loss')
    ax.grid(True)

    # ----------------------------------------------------------------------------------------------------------------------
    # Values - best set
    # ----------------------------------------------------------------------------------------------------------------------
    ax = axs[1, k]
    ax.set_title('Best set vals')

    case_dict = filtered_oparams_df  # pd.read_csv(param_file)
    case_dict['RNAP_dist'] = case_dict['RNAP_dist'].apply(lambda x: x / 1000)
    case_dict['fold_change'] = case_dict['fold_change'].apply(lambda x: x / 100)

    # Group by cluster and calculate averages and standard deviations
    cluster_stats = case_dict.groupby('cluster').agg(['mean', 'std'])

    # Extract means and standard deviations for plotting
    cluster_means = cluster_stats.loc[:, (slice(None), 'mean')].droplevel(1, axis=1)
    cluster_stds = cluster_stats.loc[:, (slice(None), 'std')].droplevel(1, axis=1)

    # Plot bar plot with error bars
    x = np.arange(len(cluster_means.columns))  # Indices for features
    width0 = 0.2  # Width of the bars

    for i, cluster in enumerate(cluster_means.index):
        ax.bar(
            x + i * width0,
            cluster_means.loc[cluster],
            width=width0,
            label=f"Cluster {cluster}",
            yerr=cluster_stds.loc[cluster],  # Add standard deviations as error bars
            capsize=5  # Add caps to the error bars
        )

    # Formatting the plot
    ax.set_title("Cluster Averages with Standard Deviations")
    ax.set_ylabel("Average Value")
    ax.set_xlabel("Features")
    ax.set_xticks(x + width0 / 2)  # Centre ticks
    ax.set_xticklabels(cluster_means.columns)
    ax.legend(title="Cluster")

    # Group by cluster and calculate averages
#    cluster_averages = case_dict.groupby('cluster').mean()
    # Plot the bar plot

    #cluster_averages.T.plot(kind='bar',ax=ax)#, palette="viridis")  # Transpose for features on the x-axis
    #ax.set_title("Cluster Averages")
    #ax.set_ylabel("Average Value")
    #ax.set_xlabel("Features")
    ##plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc='upper left')  # Adjust legend position
    ##plt.tight_layout()
plt.show()
