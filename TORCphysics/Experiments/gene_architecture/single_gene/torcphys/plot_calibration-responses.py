import pandas as pd
import matplotlib.pyplot as plt
from TORCphysics import binding_model as bm
from TORCphysics import utils
import numpy as np
# Description
#-----------------------------------------------------------------------------------------------------------------------
# The idea is to calibrate the responses of each promoter resulted from the calibration (random search)

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
dt=1.0#0.5

#k_weak=0.02
k_weak=0.0334

n_inner_workers = 9
n_subsets = 2
n_sims = n_inner_workers * n_subsets  # This is the number of simulations launched to calculate unbinding rates

model_code = 'GB-Stages-'

experimental_files = []
calibration_files = []
for pcase in promoter_cases:

    calibration_files.append('calibrate_inferred-rates/'+model_code+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.csv')

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
titles = ['weak', 'medium', 'strong']

colors = ['green', 'blue', 'red']

# Processing functions
#-----------------------------------------------------------------------------------------------------------------------
def open_rate(superhelical, k_open, threshold, width):
    U = utils.opening_function(superhelical, threshold, width)
    rate = k_open * np.exp(-U)
    return rate
# Load
#-----------------------------------------------------------------------------------------------------------------------
responses = []
for file in calibration_files:
    data = pd.read_csv(file)
    data = data.to_dict(orient='records')[0]  # Get the first (and only) dictionary
    responses.append(data)

# Plot
#-----------------------------------------------------------------------------------------------------------------------
sigma = np.arange(-.15, .0, 0.001)
# Let's plot as we do the process
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True)#, sharex=True)

x = sigma
for i, resp_dict in enumerate(responses):

    # axs[i].set_title(titles[i])
    file = calibration_files[i]
    print(i,file)

    # Gaussian binding
    gy = bm.GaussianBinding(**resp_dict).rate_modulation(superhelical=x)
    axs[0].plot(x, gy, lw=lw, color=colors[i], label=promoter_cases[i])

    # Openning rate
    oy = open_rate(x, resp_dict['k_open'], resp_dict['threshold'], resp_dict['width'])
    axs[1].plot(x, oy, lw=lw, color=colors[i])#, label=promoter_cases[i])


axs[0].grid(True)
axs[1].grid(True)
axs[0].legend(loc='best')
axs[1].set_xlabel('Superhelical density')
axs[0].set_ylabel(r'Rate ($s^{-1}$)')
axs[1].set_ylabel(r'Rate ($s^{-1}$)')
axs[0].set_title('Binding rate (Gaussian modulation)')
axs[1].set_title('Open-complex formation rate (SIST)')

# Bars of values
ax = axs[2]

nokeys = ['superhelical_op', 'spread', 'width', 'threshold', 'gamma', 'kappa', 'velocity', 'stall_torque']
# Remove the keys from all dictionaries in the list
for case_dict in responses:
    for key in nokeys:
        if key in case_dict:
            del case_dict[key]

keys = list(responses[0].keys())

# Number of cases and keys
n_cases = len(responses)
n_keys = len(keys)

# Set up the positions for each bar group (x-axis)
x = np.arange(n_keys)  # Position of each group on the x-axis
width = 0.8 / n_cases  # Dynamically calculate the width of each bar


# Plot each case in the list
for i, case_dict in enumerate(responses):

    values = [case_dict[key] for key in keys]  # Get values for the current case
    ax.bar(x + i * width - width * (n_cases - 1) / 2, values, width, label=promoter_cases[i], color=colors[i])

# Add labels, title, and custom ticks
ax.set_xlabel('Rate type')
ax.set_ylabel(r'Rate ($s^{-1}$)')
ax.set_title('Transition rates')
ax.set_xticks(x)
ax.set_xticklabels(keys)
#ax.legend()
# plt.savefig(model_code+'kw'+str(k_weak)+'_dt'+str(dt)+'.png')

#plt.savefig('calibration-response-'+model_code+'kw'+str(k_weak)+'_dt'+str(dt)+'.pdf')
#plt.savefig('calibration-response-'+model_code+'kw'+str(k_weak)+'_dt'+str(dt)+'.png')
plt.show()



