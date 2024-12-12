import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns


# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
dt=1.0#0.5


model_code = 'GB-Stages-'
model_code='sus_GB-Stages-avgx2-02-'


calibration_files = []
for pcase in promoter_cases:
    # calibration_files.append('susceptibility/'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    calibration_files.append('susceptibility/reproduce-'+model_code+pcase+'_dt'+str(dt)+'_RNAPTracking_off.pkl')



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
def prepare_local_supercoiling_distribution(results_list):
    x = results_list['distance']
    for i in range(len(x)):
        sigma = results_list['local_superhelical'][i].flatten()
        n = len(sigma)
        df = pd.DataFrame({
            'distance': np.repeat(x[i],n), 'values': sigma})
        if i == 0:
            sigma_df = df
        else:
            sigma_df = pd.concat([sigma_df, df])#, axis=0)
    return sigma_df

def calculate_local_supercoiling_accumulation(results_list):
    x = results_list['distance']
    for i in range(len(x)):
        sigma = results_list['local_superhelical_variations'][i].flatten()
        n = len(sigma)
        df = pd.DataFrame({
            'distance': np.repeat(x[i],n), 'values': sigma})
        if i == 0:
            sigma_df = df
        else:
            sigma_df = pd.concat([sigma_df, df])
    return sigma_df

def calculate_global_supercoiling_accumulation(results_list):
    x = results_list['distance']
    for i in range(len(x)):
        sigma = results_list['global_superhelical_variations'][i].flatten()
        n = len(sigma)
        df = pd.DataFrame({
            'distance': np.repeat(x[i],n), 'values':sigma})
        if i == 0:
            sigma_df = df
        else:
            sigma_df = pd.concat([sigma_df, df])
    return sigma_df

def prepare_global_supercoiling_distribution(results_list):
    x = results_list['distance']
    for i in range(len(x)):
        sigma = results_list['global_superhelical'][i].flatten()
        n = len(sigma)
        df = pd.DataFrame({
            'distance': np.repeat(x[i],n), 'values':sigma})
        if i == 0:
            sigma_df = df
        else:
            sigma_df = pd.concat([sigma_df, df])
    return sigma_df

def prepare_active_supercoiling_distribution(results_list):
    x = results_list['distance']
    for i in range(len(x)):
        sigma = results_list['active_superhelical'][i].flatten()
        n = len(sigma)
        df = pd.DataFrame({
            'distance': np.repeat(x[i],n), 'values':sigma})
        if i == 0:
            sigma_df = df
        else:
            sigma_df = pd.concat([sigma_df, df])
    return sigma_df


# Load
#-----------------------------------------------------------------------------------------------------------------------
pickle_data = []
for pickle_file in calibration_files:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_data.append(data)

active_sigma = []
for i, data in enumerate(pickle_data):
    x=i
    active_sigma.append(prepare_active_supercoiling_distribution(data[0]))

local_sigma = []
for i, data in enumerate(pickle_data):
    x=i
    local_sigma.append(prepare_local_supercoiling_distribution(data[0]))

global_sigma = []
for i, data in enumerate(pickle_data):
    x=i
    global_sigma.append(prepare_global_supercoiling_distribution(data[0]))

local_dsigma = []
for i, data in enumerate(pickle_data):
    x=i
    local_dsigma.append(calculate_local_supercoiling_accumulation(data[0]))#

global_dsigma = []
for i, data in enumerate(pickle_data):
    x=i
    global_dsigma.append(calculate_global_supercoiling_accumulation(data[0]))

# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we do the process
fig, axs = plt.subplots(3,2, figsize=(2*width, 3*height), tight_layout=True, sharex=False)

# Local
var = local_sigma
#var = local_dsigma
for i, sigma in enumerate(var):
    ax = axs[i,0]

    ax.set_title(titles[i], color=colors[i])

    # Check for duplicate values in the 'distance' column
    # sigma = sigma[sigma.duplicated(subset=['distance'], keep=False)]
    # Filter for duplicate 'distance' values
    sigma = sigma[sigma.duplicated(subset=['distance'], keep=False)]

    # Reset the index to ensure uniqueness
    sigma = sigma.reset_index(drop=True)

    sns.violinplot(x='distance', y='values', data=sigma, ax=ax, inner="quart", cut=0, color=colors[i])

    ax.set_ylabel(r'Local superhelicity')
    ax.set_xlabel(r'Upstream distance')
    ax.set_ylim([-.1, -0.0])
    #ax.set_ylim([-0.005, 0.005])
    #ax.grid(True)


# Global
var = global_sigma
#var = global_dsigma
for i, sigma in enumerate(var):
    ax = axs[i,1]

    ax.set_title(titles[i], color=colors[i])

    # Check for duplicate values in the 'distance' column
    sigma = sigma[sigma.duplicated(subset=['distance'], keep=False)]

    # Reset the index to ensure uniqueness
    sigma = sigma.reset_index(drop=True)

    sns.violinplot(x='distance', y='values', data=sigma, ax=ax, inner="quart", cut=0, color=colors[i])

    ax.set_ylabel(r'Global superhelicity')
    ax.set_xlabel(r'Upstream distance')
    ax.set_ylim([-.1, -0.0])
    # ax.set_ylim([-0.0005, 0.0005])
    #ax.grid(True, zorder=10)

plt.show()



outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

# Let's plot as we do the process
fig, axs = plt.subplots(3,2, figsize=(2*width, 3*height), tight_layout=True, sharex=False)

# Local
#var = local_sigma
var = local_dsigma
#var = global_sigma
for i, sigma in enumerate(var):
    ax = axs[i,0]

    ax.set_title(titles[i], color=colors[i], fontsize=title_size)

    # Check for duplicate values in the 'distance' column
    sigma = sigma[sigma.duplicated(subset=['distance'], keep=False)]

    # Reset the index to ensure uniqueness
    sigma = sigma.reset_index(drop=True)

    sns.violinplot(x='distance', y='values', data=sigma, ax=ax, inner="quart", cut=0, color=colors[i])

    # Add label outside the plot
    ax.text(-0.1, 1.1, outside_label[i], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

    ax.set_ylabel(r'Global superhelicity', fontsize=font_size)
    ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    ax.set_ylim([-.1, -0.0])
    #ax.set_ylim([-0.005, 0.005])
    ax.grid(True)


# Global
#var = global_sigma
# var = global_dsigma
var = active_sigma
for i, sigma in enumerate(var):
    ax = axs[i,1]

    ax.set_title(titles[i], color=colors[i], fontsize=title_size)

    # Check for duplicate values in the 'distance' column
    sigma = sigma[sigma.duplicated(subset=['distance'], keep=False)]

    # Reset the index to ensure uniqueness
    sigma = sigma.reset_index(drop=True)

    sns.violinplot(x='distance', y='values', data=sigma, ax=ax, inner="quart", cut=0, color=colors[i])

    # Add label outside the plot
    ax.text(-0.1, 1.1, outside_label[i+3], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

    ax.set_ylabel(r'Active superhelicity', fontsize=font_size)
    ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    ax.set_ylim([-.1, 0.0])
    #ax.set_ylim([-0.0005, 0.0005])
    ax.grid(True, zorder=10)

#plt.savefig('supercoiling_'+model_code+'kw'+str(k_weak)+'_dt'+str(dt)+'.png')
#plt.savefig('supercoiling_'+model_code+'kw'+str(k_weak)+'_dt'+str(dt)+'.pdf')

plt.show()
