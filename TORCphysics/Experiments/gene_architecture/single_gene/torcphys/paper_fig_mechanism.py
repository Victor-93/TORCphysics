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
outfile = model_code+'mechanism'

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

    #dist_files_tracking_on.append('susceptibility/dist_reproduce-'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    dist_files_tracking_on.append('susceptibility/dist_reproduce-'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    dist_files_tracking_off.append('susceptibility/dist_reproduce-'+model_code+pcase+'_dt'+str(dt)+'_RNAPTracking_off.pkl')

dist_files_tracking_on = calibration_files_tracking_on
dist_files_tracking_off = calibration_files_tracking_off
# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 7
height = 3.75
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
def prepare_dists(x, dist_list, label):
    combined = []
    for distance, df in zip(x, dist_list):
        # Melt the DataFrame to long format
        melted = df.melt(var_name="simulation", value_name="value")
        melted["distance"] = distance  # Add the distance column
        melted["type"] = label  # Add a label for on/off
        combined.append(melted)
    return pd.concat(combined, ignore_index=True)
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
fig, axs = plt.subplots(2,2, figsize=(2*width, 2*height), tight_layout=True, sharex=True)

outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

color_on = 'black'
color_off = 'red'

color_enzyme = {'RNAP': 'black', 'topoI':'red', 'gyrase':'blue'}
nenzyme_ylim = [0,6]

n = len(calibration_files_tracking_on)
for i in range(n):

    # Susceptibility
    # ----------------------------------------------------------------------------------
    ax = axs[0,0]
    array_on = susceptibility_tracking_on[i]
    array_off = susceptibility_tracking_off[i]

    ax.set_title('Susceptibility ', fontsize=title_size)

    # Get variables
    # Tracking ON
    x_on = array_on[0]
    y_on = array_on[1]
    ys_on = array_on[2]

    # Tracking OFF
    x_off = array_off[0]
    y_off = array_off[1]
    ys_off = array_off[2]

    ax.plot(x_on, y_on, model_ls, lw=lw, color=color_on, label='ON')
    ax.plot(x_off, y_off, model_ls, lw=lw, color=color_off, label='OFF')

    ax.fill_between(x_on, y_on-ys_on, y_on+ys_on, lw=lw, color=color_on, alpha=0.2)
    ax.fill_between(x_off, y_off-ys_off, y_off+ys_off, lw=lw, color=color_off, alpha=0.2)

    #ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    ax.set_ylabel(r'Susceptiblity', fontsize=font_size)
    ax.grid(True)


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

    ax.plot(x_on, y_on, model_ls, lw=lw, color=color_on, label='ON')
    ax.plot(x_off, y_off, model_ls, lw=lw, color=color_off, label='OFF')

    ax.fill_between(x_on, y_on-ys_on, y_on+ys_on, lw=lw, color=color_on, alpha=0.2)
    ax.fill_between(x_off, y_off-ys_off, y_off+ys_off, lw=lw, color=color_off, alpha=0.2)

    ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    ax.set_ylabel(r'Production rate $(s^{-1})$', fontsize=font_size)
    ax.grid(True)
    ax.set_xlim([80, 6000])
    ax.set_xscale('log')

    # THIRD OPTION: Plot distributions
    # SECOND OPTION: Plot RNAPs and topoisomerases together
    # Number of bound Topoisomerase I
    # ----------------------------------------------------------------------------------
    ax = axs[0,1]
    data_on = pickle_dist_tracking_on[i][0]
    data_off = pickle_dist_tracking_off[i][0]
    name = 'topoI'

    ax.set_title('Bound Topoisomerase I', fontsize=title_size) #color=colors[i],

    # Get data
    x = np.array(data_on['distance'])
    y_on= []
    ys_on = []
    for j in range(len(x)):
        y_on.append(np.mean([d for d in data_on['nenzymes'][name][j]]))
        ys_on.append(np.std([d for d in data_on['nenzymes'][name][j]]))

    y_off= []
    ys_off = []
    for j in range(len(x)):
        y_off.append(np.mean([d for d in data_off['nenzymes'][name][j]]))
        ys_off.append(np.std([d for d in data_off['nenzymes'][name][j]]))

    ax.errorbar(x,y_on,yerr=ys_on, fmt='-o', color=color_on, label=name+' ON')
    ax.errorbar(x,y_off,yerr=ys_off, fmt='-o', color=color_off, label=name+' OFF')
    #ax.plot(x, y_on, model_ls, lw=lw, color=color_on, label=name + ' ON')
    #ax.plot(x, y_off, model_ls, lw=lw, color=color_off, label=name + ' OFF')

    # Gyrase additional
    plot_gyrase = True
    if plot_gyrase:
        name = 'gyrase'

        y_on = []
        ys_on = []
        for j in range(len(x)):
            y_on.append(np.mean([d for d in data_on['nenzymes'][name][j]]))
            ys_on.append(np.std([d for d in data_on['nenzymes'][name][j]]))

        y_off = []
        ys_off = []
        for j in range(len(x)):
            y_off.append(np.mean([d for d in data_off['nenzymes'][name][j]]))
            ys_off.append(np.std([d for d in data_off['nenzymes'][name][j]]))

        ax.errorbar(x, y_on, yerr=ys_on, fmt='-o', color='blue', label='gyrase ON')
        ax.errorbar(x, y_off, yerr=ys_off, fmt='-o', color='cyan', label='gyrase OFF')
        #ax.plot(x,y_on, model_ls, lw=lw,color='blue', label=name+' ON')
        #ax.plot(x,y_off, model_ls, lw=lw,color='cyan', label=name+' OFF')

    ax.set_ylabel(r'Number of bound Topoisomerases', fontsize=font_size)
    ax.grid(True)
    #ax.set_ylim(nenzyme_ylim)

    # Number of bound RNAPs
    # ----------------------------------------------------------------------------------
    ax = axs[1,1]
    name = 'RNAP'

    ax.set_title('Bound RNAPs', fontsize=title_size) #color=colors[i],

    # Get data
    x = np.array(data_on['distance'])
    y_on= []
    ys_on = []
    for j in range(len(x)):
        y_on.append(np.mean([d[:2500] for d in data_on['nenzymes'][name][j]]))
        ys_on.append(np.std([d[:2500] for d in data_on['nenzymes'][name][j]]))

    y_off= []
    ys_off = []
    for j in range(len(x)):
        y_off.append(np.mean([d[:2500] for d in data_off['nenzymes'][name][j]]))
        ys_off.append(np.std([d[:2500] for d in data_off['nenzymes'][name][j]]))

    ax.errorbar(x,y_on,yerr=ys_on, fmt='-o', color=color_on, label=name+' ON')
    ax.errorbar(x,y_off,yerr=ys_off, fmt='-o', color=color_off, label=name+' OFF')

    #ax.plot(x,y_on, model_ls, lw=lw,color=color_on, label=name+' ON')
    #ax.plot(x,y_off, model_ls, lw=lw,color=color_off, label=name+' OFF')

    ax.set_ylabel(r'Number of bound RNAPs', fontsize=font_size)
    ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    ax.grid(True)
    ax.set_xlim([80, 6000])
    ax.set_xscale('log')


ax_list = [axs[0,0], axs[1,0], axs[0,1], axs[1,1]]
for i, ax in enumerate(ax_list):
    ax.legend(loc='upper left')

    # Add label outside the plot
    ax.text(-0.12, 1., outside_label[i], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

plt.savefig(outfile+'.pdf')
plt.show()
