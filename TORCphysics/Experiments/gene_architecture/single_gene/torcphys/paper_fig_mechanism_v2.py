import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
from TORCphysics import utils
from TORCphysics import binding_model as bm

# Description
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot susceptibility/production rate in the same panel (2 axes) when the tracking is on/off, and
# maybe include the distribution of bound enzymes of topoI and RNAP

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
promoter_cases = ['medium']
dt=1.0#0.5

model_code='sus_GB-Stages-avgx2-02-'
outfile = model_code+'fig'

susceptibility_files = []
calibration_files_tracking_on = []  # Topo I follows RNAPs
calibration_files_tracking_off = []  # Topo I does not follow RNAPs
dist_files_tracking_on =[] # This one have the distributions
dist_files_tracking_off =[]
params_files = []
for pcase in promoter_cases:
    susceptibility_files.append('../junier_data/'+ pcase + '.csv')

    calibration_files_tracking_on.append('susceptibility/'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    calibration_files_tracking_off.append('susceptibility/reproduce-'+model_code+pcase+'_dt'+str(dt)+'_RNAPTracking_off.pkl')
    params_files.append('table_'+model_code+pcase+'_dt'+str(dt)+'.csv')

    dist_files_tracking_on.append('susceptibility/dist_reproduce-'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    dist_files_tracking_off.append('susceptibility/dist_reproduce-'+model_code+pcase+'_dt'+str(dt)+'_RNAPTracking_off.pkl')

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 7
height = 3.5
lw = 4
font_size = 12
xlabel_size = 14
title_size = 16

# line styles
model_ls = '-o'
exp_ls = '--o'
titles = ['Weak Promoter', 'Medium Promoter', 'Strong Promoter']

colors = ['green', 'blue', 'red']

# Processing functions
#-----------------------------------------------------------------------------------------------------------------------
# Calculate rates as a function of distance
def get_prod_rates(results_list):
    x = results_list['distance'] #[d['distance'] for d in results_list]
    y = []
    ys = []
    # List with just results (not distance and not dict anymore)
    for j, rates_array in enumerate(results_list['prod_rate']): #case_results is the rate

        # Calculate mean and standard deviation
        mean = rates_array[0]
        std = rates_array[1]
        y.append(mean)
        ys.append(std)

    rates = np.array([x, y, ys])
    return rates

# Calculate susceptibilities as a function of distance
def get_susceptibility(results_list):
    x = results_list['distance'] #[d['distance'] for d in results_list]
    y = []
    ys = []
    # List with just results (not distance and not dict anymore)
    for j, sus_array in enumerate(results_list['susceptibility']): #case_results is the rate

        # Calculate mean and standard deviation
        mean = sus_array[0]
        std = sus_array[1]
        y.append(mean)
        ys.append(std)

    rates = np.array([x, y, ys])
    return rates


def open_rate(superhelical, k_open, threshold, width):
    U = utils.opening_function(superhelical, threshold, width)
    rate = k_open * np.exp(-U)
    return rate

# Prepare data for plotting
# data = All pkl data. label (ON or OFF), ename=enzyme name e.g., topoI
def prepare_dists(data, label, ename):
    x = data['distance']
    n = len(x)
    dist_list = []
    for k in range(n):
        dist_list.append(data['nenzymes_dist'][ename][k])

    combined = []
    for distance, df in zip(x, dist_list):
        # Melt the DataFrame to long format
        melted = df.melt(var_name="simulation", value_name="value")
        melted["distance"] = distance  # Add the distance column
        melted["type"] = label  # Add a label for on/off
        combined.append(melted)
    output = pd.concat(combined, ignore_index=True)
    return output

# Load
#-----------------------------------------------------------------------------------------------------------------------
pickle_data_tracking_on = []
for pickle_file in calibration_files_tracking_on:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_data_tracking_on.append(data)

pickle_data_tracking_off = []
for pickle_file in calibration_files_tracking_off:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_data_tracking_off.append(data)

pickle_dist_tracking_on = []
for pickle_file in dist_files_tracking_on:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_dist_tracking_on.append(data)

pickle_dist_tracking_off = []
for pickle_file in dist_files_tracking_off:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_dist_tracking_off.append(data)

# Calculate rates
rates_tracking_on = []
for i, data in enumerate(pickle_data_tracking_on):
    x=i
    rates_tracking_on.append(get_prod_rates(data[0]))

# Calculate rates
rates_tracking_off = []
for i, data in enumerate(pickle_data_tracking_off):
    x=i
    rates_tracking_off.append(get_prod_rates(data[0]))

# Calculate susceptibilitie
susceptibility_tracking_on = []
for i, data in enumerate(pickle_data_tracking_on):
    x=i
    susceptibility_tracking_on.append(get_susceptibility(data[0]))

# Calculate susceptibilitie
susceptibility_tracking_off = []
for i, data in enumerate(pickle_data_tracking_off):
    x=i
    susceptibility_tracking_off.append(get_susceptibility(data[0]))

responses = []
for file in params_files:
    data = pd.read_csv(file)
    data = data.to_dict(orient='records')[0]  # Get the first (and only) dictionary
    responses.append(data)

# Plot
#-----------------------------------------------------------------------------------------------------------------------
sigma = np.arange(-.13, .0, 0.001)
#sigma = np.arange(-.2, .0, 0.001)
# Let's plot as we do the process
fig, axs = plt.subplots(2,2, figsize=(2*width, 2*height), tight_layout=True)#, sharex=True)

outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

color_on = 'black'
color_off = 'green'

color_enzyme = {'RNAP': 'black', 'topoI':'red', 'gyrase':'blue'}
nenzyme_ylim = [0,6]

n = len(calibration_files_tracking_on)
for i in range(n):

    # Susceptibility
    # ----------------------------------------------------------------------------------
    ax = axs[0,0]
    array_on = susceptibility_tracking_on[i]
    array_off = susceptibility_tracking_off[i]

    ax.set_title('Susceptibility', fontsize=title_size)

    # Get variables
    # Tracking ON
    x_on = array_on[0]
    y_on = array_on[1]
    ys_on = array_on[2]

    # Tracking OFF
    x_off = array_off[0]
    y_off = array_off[1]
    ys_off = array_off[2]

    ax.plot(x_on, y_on, model_ls, lw=lw, color=color_on, label='On')
    ax.plot(x_off, y_off, model_ls, lw=lw, color=color_off, label='On')

    ax.fill_between(x_on, y_on-ys_on, y_on+ys_on, lw=lw, color=color_on, alpha=0.2)
    ax.fill_between(x_off, y_off-ys_off, y_off+ys_off, lw=lw, color=color_off, alpha=0.2)

    #ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    ax.set_ylabel(r'Susceptiblity', fontsize=font_size)
    ax.grid(True)
    ax.set_xscale('log')


    # Production rate
    # ----------------------------------------------------------------------------------
    ax = axs[1,0]
    array_on = rates_tracking_on[i]
    array_off = rates_tracking_off[i]

    ax.set_title('Production rate', fontsize=title_size) #color=colors[i],

    # Get variables
    # Tracking ON
    x_on = array_on[0]
    y_on = array_on[1]
    ys_on = array_on[2]

    # Tracking OFF
    x_off = array_off[0]
    y_off = array_off[1]
    ys_off = array_off[2]

    ax.plot(x_on, y_on, model_ls, lw=lw, color=color_on, label='On')
    ax.plot(x_off, y_off, model_ls, lw=lw, color=color_off, label='On')

    ax.fill_between(x_on, y_on-ys_on, y_on+ys_on, lw=lw, color=color_on, alpha=0.2)
    ax.fill_between(x_off, y_off-ys_off, y_off+ys_off, lw=lw, color=color_off, alpha=0.2)

    ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    ax.set_ylabel(r'Production rate $(s^{-1})$', fontsize=font_size)
    ax.grid(True)
    ax.set_xlim([80, 6000])
    ax.set_xscale('log')

    # THIRD OPTION: Plot distributions
    # Number of bound Topoisomerase I
    # ----------------------------------------------------------------------------------
    ax = axs[0,1]
    data_on = pickle_dist_tracking_on[i][0]
    data_off = pickle_dist_tracking_off[i][0]
    name = 'topoI'

    ax.set_title('Bound Topoisomerase I', fontsize=title_size) #color=colors[i],

    # Combine data for ON and OFF
    data_on = prepare_dists(data_on, "ON", name)
    data_off = prepare_dists(data_off, "OFF", name)
    all_data = pd.concat([data_on, data_off], ignore_index=True)

    # Plot distributions as a function of distance
    #plt.figure(figsize=(12, 6))
    sns.violinplot(x="distance", y="value", hue="type", data=all_data, split=True, inner="quart", ax=ax)
#    plt.title("Distributions of Values as a Function of Distance")
    ax.legend(title="Type")

    ax.set_ylabel(r'Number of bound Topoisomerase I', fontsize=font_size)
    ax.grid(True)
    ax.set_ylim([-1,10])

    # Number of bound RNAP
    # ----------------------------------------------------------------------------------
    ax = axs[1,1]
    data_on = pickle_dist_tracking_on[i][0]
    data_off = pickle_dist_tracking_off[i][0]
    name = 'RNAP'

    ax.set_title('Bound RNAP', fontsize=title_size) #color=colors[i],

    # Combine data for ON and OFF
    data_on = prepare_dists(data_on, "ON", name)
    data_off = prepare_dists(data_off, "OFF", name)
    all_data = pd.concat([data_on, data_off], ignore_index=True)

    # Plot distributions as a function of distance
    #plt.figure(figsize=(12, 6))
    sns.violinplot(x="distance", y="value", hue="type", data=all_data, split=True, inner="quart", ax=ax)
#    plt.title("Distributions of Values as a Function of Distance")
#    ax.legend(title="Type")

    ax.set_ylabel(r'Number of bound RNAPs', fontsize=font_size)
    ax.grid(True)
    ax.set_ylim([-1,4])

plt.show()

