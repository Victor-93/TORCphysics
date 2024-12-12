import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns


# Let's just plot supercoilng as error bars
#-----------------------------------------------------------------------------------------------------------------------


# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['medium']
dt=1.0#0.5


model_code = 'GB-Stages-'
model_code='sus_GB-Stages-avgx2-02-'


calibration_files_tracking_on = []  # Topo I follows RNAPs
calibration_files_tracking_off = []  # Topo I does not follow RNAPs
for pcase in promoter_cases:
    calibration_files_tracking_on.append('susceptibility/'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    calibration_files_tracking_off.append('susceptibility/reproduce-'+model_code+pcase+'_dt'+str(dt)+'_RNAPTracking_off.pkl')

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 7
height = 3
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

# line styles
model_ls = '-o'
exp_ls = '--o'
titles = ['Weak promoter', 'Medium promoter', 'Strong promoter']

colors = ['green', 'blue', 'red']

# Processing functions
#-----------------------------------------------------------------------------------------------------------------------
def get_sigma(results_list, sigma_type):
    x = results_list['distance']
    y = [] # average
    ys = [] #standard deviation
    for i in range(len(x)):
        sigma = results_list[sigma_type][i].flatten()
        y.append( np.mean(sigma) )
        ys.append( np.std(sigma) )
    x = np.array(x)
    y = np.array(y)
    ys = np.array(ys)
    return x,y,ys

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
# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we do the process
#fig, axs = plt.subplots(3,figsize=(width, 3*height), tight_layout=True, sharex=False)
fig, axs = plt.subplots(2,figsize=(width, 2*height), tight_layout=True, sharex=False)

n = len(calibration_files_tracking_on)
color_on = 'black'
color_off = 'red'

variables = ['local_superhelical', 'global_superhelical']#, 'active_superhelical']
title_list = ['Local superhelicity', 'Global superhelicity', 'Active superhelicity']
ylabel_list = ['Local Superhelicity', 'Global Superhelicity', 'Active superhelicity']
for i, var in enumerate(variables):
    # Local superhelical
    ax = axs[i]
    #ax.set_title(title_list[i], fontsize=title_size)

    data_on = pickle_data_tracking_on[0][0]
    data_off = pickle_data_tracking_off[0][0]

    x, y_on, ys_on = get_sigma(data_on, var)
    x, y_off, ys_off = get_sigma(data_off, var)

    ax.errorbar(x,y_on,yerr=ys_on, fmt='-o', color=color_on, label='ON')
    ax.errorbar(x,y_off,yerr=ys_off, fmt='-o', color=color_off, label='OFF')

    ax.set_ylabel(ylabel_list[i], fontsize=font_size)
    ax.grid(True)
    ax.set_xlim([80, 6000])
    ax.set_xscale('log')

axs[1].set_xlabel('Upstream Distance (bp)', fontsize=font_size)

plt.show()
