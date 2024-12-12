import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns


# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
dt=1.0#0.5

k_weak=0.025
# k_weak=0.0334

model_code = 'GB-Stages-'
model_code='sus_GB-Stages-avgx2-02-'


experimental_files = []
calibration_files = []
for pcase in promoter_cases:
    experimental_files.append('../junier_data/inferred-rate_kw'+str(k_weak)+'_' + pcase + '.csv')
    calibration_files.append('calibrate_inferred-rates/reproduce-'+model_code+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.pkl')
    calibration_files.append('susceptibility/'+model_code+pcase+'_dt'+str(dt)+'.pkl')


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
    x = [d['distance'] for d in results_list]
    # y = np.zeros_like(x)
    # ys = np.zeros_like(x)
    local_sigma = [d['result']['local_superhelical'] for d in results_list]

    n = len(local_sigma[0][0])

    for i in range(len(x)):
        df = pd.DataFrame({
            'distance': np.repeat(x[i],n), 'values': local_sigma[i][0]})
        if i == 0:
            sigma_df = df
        else:
            sigma_df = pd.concat([sigma_df, df])#, axis=0)

    return sigma_df

def calculate_local_supercoiling_accumulation(results_list):
    x = [d['distance'] for d in results_list]

    sigma_list = [d['result']['local_superhelical'] for d in results_list]
    n = len(sigma_list[0][0])

    for j, sigma_array in enumerate(sigma_list): #case_results is the rate
        dsigma = np.zeros(n-1)
        for i in range(n-1):
            dsigma[i] = (sigma_array[0][i+1] - sigma_array[0][i])/dt

        df = pd.DataFrame({
            'distance': np.repeat(x[j], n-1), 'values': dsigma})
        if j == 0:
            dsigma_df = df
        else:
            dsigma_df = pd.concat([dsigma_df, df])  # , axis=0)

    return dsigma_df

def calculate_global_supercoiling_accumulation(results_list):
    x = [d['distance'] for d in results_list]

    sigma_list = [d['result']['global_superhelical'] for d in results_list]
    n = len(sigma_list[0][0])

    for j, sigma_array in enumerate(sigma_list): #case_results is the rate
        dsigma = np.zeros(n-1)
        for i in range(n-1):
            dsigma[i] = (sigma_array[0][i+1] - sigma_array[0][i])/dt

        df = pd.DataFrame({
            'distance': np.repeat(x[j], n-1), 'values': dsigma})
        if j == 0:
            dsigma_df = df
        else:
            dsigma_df = pd.concat([dsigma_df, df])  # , axis=0)

    return dsigma_df

def prepare_global_supercoiling_distribution(results_list):
    x = [d['distance'] for d in results_list]
    # y = np.zeros_like(x)
    # ys = np.zeros_like(x)
    global_sigma = [d['result']['global_superhelical'] for d in results_list]

    n = len(global_sigma[0][0])

    for i in range(len(x)):
        df = pd.DataFrame({
            'distance': np.repeat(x[i],n), 'values': global_sigma[i][0]})
        if i == 0:
            sigma_df = df
        else:
            sigma_df = pd.concat([sigma_df, df])#, axis=0)

    return sigma_df

# Load
#-----------------------------------------------------------------------------------------------------------------------
pickle_data = []
for pickle_file in calibration_files:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_data.append(data)

local_sigma = []
for i, data in enumerate(pickle_data):
    x=i
    local_sigma.append(prepare_local_supercoiling_distribution(data[0]['data']))

global_sigma = []
for i, data in enumerate(pickle_data):
    x=i
    global_sigma.append(prepare_global_supercoiling_distribution(data[0]['data']))

local_dsigma = []
for i, data in enumerate(pickle_data):
    x=i
    local_dsigma.append(calculate_local_supercoiling_accumulation(data[0]['data']))

global_dsigma = []
for i, data in enumerate(pickle_data):
    x=i
    global_dsigma.append(calculate_global_supercoiling_accumulation(data[0]['data']))

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
    sigma = sigma[sigma.duplicated(subset=['distance'], keep=False)]

    sns.violinplot(x='distance', y='values', data=sigma, ax=ax, inner="quart", cut=0, color=colors[i])

    ax.set_ylabel(r'Local superhelicity')
    ax.set_xlabel(r'Upstream distance')
    ax.set_ylim([-.06, -0.03])
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

    sns.violinplot(x='distance', y='values', data=sigma, ax=ax, inner="quart", cut=0, color=colors[i])

    ax.set_ylabel(r'Global superhelicity')
    ax.set_xlabel(r'Upstream distance')
    ax.set_ylim([-.06, -0.03])
    # ax.set_ylim([-0.0005, 0.0005])
    #ax.grid(True, zorder=10)

plt.show()



