import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns
from TORCphysics import utils
from TORCphysics import binding_model as bm
import sys

# Description
#-----------------------------------------------------------------------------------------------------------------------
# Lets try to define the new loss, by analysing the function and its components.

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
#promoter_cases = ['medium']
promoter_cases = ['weak', 'medium', 'strong']
promoter_labels = ['Weak', 'Medium', 'Strong']
#variables = ['loss', 'prod_rate', 'local_superhelical', 'global_superhelical']
#colors = ['green', 'blue', 'red']
percentage_threshold = 0.05
outfig = 'genearch_responses'
mean_weight = 1.0
error_weight = 1.0#1#0.1#1.0


# This list of dicts will include all the input information
info_list = [
    # genearch_V0 - gamma 0.05
    {**{pm: f"../calibration/sus_genearch_V0-01-{pm}_gamma0.05_dt1.0-trials.pkl" for pm in promoter_cases},
     'label': 'V0'},

    # genearch_V1 - gamma 0.1
    {**{pm: f"../calibration/sus_genearch_V1-01-{pm}_gamma0.1_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V1'},

    # genearch_V2
    {**{pm: f"../calibration/sus_genearch_V2-01-{pm}_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V2'}

]

susceptibility_files = []
promoter_responses_files = []
for pcase in promoter_cases:
    susceptibility_files.append('../junier_data/'+ pcase + '.csv')

    # This contains the promoter responses
    promoter_responses_files.append('../promoter_responses/' + pcase + '.csv')

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 5
height = 3#3.75
lw = 4
font_size = 12
xlabel_size = 14
title_size = 16
ms=6

# line styles
model_ls = '-o'
exp_ls = '--o'

colors = ['green', 'blue', 'red', 'purple', 'orange']

# Functions
#-----------------------------------------------------------------------------------------------------------------------
def calculate_susceptibility(data_trials):
    for trial in data_trials:
        pr = trial['prod_rate'][0]
        prs = trial['prod_rate'][1]

        y = pr/pr[4] #susceptibility
        ys = np.zeros_like(y)
        for j in range(len(y)):
            ys[j] = y[j] * np.sqrt( np.square(prs[j]/pr[j]) + np.square(prs[4]/pr[4]) )

    return y, ys

def get_loss(data_trials, exp_mean,exp_err, w_mean=1, w_err=1):
    # w is the weight for the mean and ws for the standard deviation

    loss_list = []
    n = len(exp_mean)
    # First calculate susceptibility
    for trial in data_trials:
        pr = trial['prod_rate'][:,0]
        prs = trial['prod_rate'][:,1]

        y = pr / pr[4]  # susceptibility
        ys = np.zeros_like(y)

        # Sometimes prj
        for j in range(len(y)):
            ys[j] = y[j] * np.sqrt(np.square(prs[j] / pr[j]) + np.square(prs[4] / pr[4]))

        loss = 0
        for i in range(n):
            loss += w_mean*np.square(exp_mean[i] - y[i]) + w_err*np.square(exp_err[i] - ys[i])
        loss_list.append(loss)
    return loss_list

def open_rate(superhelical, k_open, threshold, width):
    U = utils.opening_function(superhelical, threshold, width)
    rate = k_open * np.exp(-U)
    return rate

# Process and plot - Let's plot susceptibility vs upstream barrier distance of best cases
#-----------------------------------------------------------------------------------------------------------------------
nokeys = ['superhelical_op', 'spread', 'width', 'threshold', 'gamma', 'kappa', 'velocity', 'stall_torque']
outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']
ncols = 2
nrows = 3
outfig = 'genearch_responses'
pelE_color = 'gray'
spacer_color = 'black'
t_rates_xlim = [0,0.51]
#1 SIST, 2 Gaussian V1, 3 Gaussian V2. 4 responses V0, 5 responses V1, 6 Responses V2
fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*width, nrows*height),
tight_layout=True)#, sharex=True, sharey=True)


sigma = np.arange(-.13, .01, 0.001)

# a) SIST Responses
# ----------------------------------------------------
ax = axs[0,0]
olabel = outside_label[0]

# Load
responses = []
for file in promoter_responses_files:
    data = pd.read_csv(file)
    data = data.to_dict(orient='records')[0]  # Get the first (and only) dictionary
    responses.append(data)

# Process and plot
x = sigma
# Now let's plot the responses
for i in range(3):

    # NOTE: I am using a k_open of 1.0 because it doesn't matter as we are going to normalize it.
    #       We don't need it's value right now.
    # Openning rate
    oy = open_rate(x, 1.0, responses[i]['threshold'], responses[i]['width'])
    y =  (oy - np.min(oy)) / (np.max(oy) - np.min(oy)) # Let's normalize
    ax.plot(x, y, lw=lw, color=colors[i], label=promoter_cases[i])#promoter_labels[i])

# Meyer model for PelE? model
mb= bm.MeyerPromoterOpening (**{'k_on':1.0})
y = mb.rate_modulation(superhelical=x)
y = (y - np.min(y)) / (np.max(y) - np.min(y))  # Let's normalize
ax.plot(x, y, lw=lw, color=pelE_color, label='pelE')

ax.text(-0.1, 1.1, olabel, transform=ax.transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

ax.grid(True)
ax.legend(loc='best',fontsize=font_size)
ax.set_xlabel('Superhelical Density', fontsize=xlabel_size)
ax.set_ylabel(r'Rate Modulation', fontsize=xlabel_size)
ax.set_title('Open Complex Formation Rate (SIST)', fontsize=title_size)

# b) Gaussian V1
# ----------------------------------------------------
V_num = 1
ax = axs[1,0]
olabel = outside_label[1]
mytitle = 'Gaussian Modulation: V1'
case = info_list[V_num]
responses_list_V1 = []
# Go through promoter
for j, pm in enumerate(promoter_cases):


    my_color = colors[j]

    figtitle = case['label'] + ': ' + promoter_labels[j] + ' Promoter'
    ax.set_title(figtitle, fontsize=title_size)  # color=colors[i],

    # Experimental susceptibility
    exp = pd.read_csv(susceptibility_files[0])  # read
    exp_x = exp['distance']
    exp_y = exp['Signal']
    exp_ys = exp['Error']

    # Load and prepare data, sort losses and find best case according to the weighted loss
    # --------------------------------------------------------------------
    with open(case[pm], 'rb') as file:
        data = pickle.load(file).results

    with open(case[pm], 'rb') as file:
        vals_dictlist = pickle.load(file).vals

    # Change std for error in the mean (or standard error)
    for dat in data:
        dat['prod_rate'][:, 1] = dat['prod_rate'][:, 1] / np.sqrt(60)

    # Filter data to include only entries with 'status' == 'ok'
    results = [item for item in data if item['status'] == 'ok']

    # Calculate new loss and just added to the list with results "results"
    for dat in results:
        new_loss = get_loss([dat], exp_y, exp_ys, w_mean=mean_weight, w_err=error_weight)
        dat['new_loss'] = new_loss[0]

    # Filter data_list to include only entries with 'status' == 'ok'
    dat_list = [item for item in results if item['status'] == 'ok']

    # Find the dictionary with the smallest 'loss'
    dat = min(dat_list, key=lambda x: x['new_loss'])

    # Get the index of this min_loss_data in the original data_list and append it
    min_loss_index = next(index for index, item in enumerate(dat_list) if item == dat)

    # Collect responses dict
    # ------------------------------------------------------------------------------
    # Let's form the response dict
    resp_dict = {}
    # Getting values from trials
    for key in vals_dictlist:
        resp_dict[key] = vals_dictlist[key][min_loss_index]

    # Getting values from promoter response
    for key in ['width', 'threshold']:
        resp_dict[key] = responses[j][key]
    # Rename kon to k_on
    resp_dict['k_on'] = resp_dict.pop('kon')

    # Let's save the responses
    responses_list_V1.append(resp_dict)

    # Plot Gaussian
    # ------------------------------------------------------------------------------
    # Gaussian binding
    gy = bm.GaussianBinding(**resp_dict).rate_modulation(superhelical=x)
    y =  (gy - np.min(gy)) / (np.max(gy) - np.min(gy)) # Let's normalize
    ax.plot(x, y, lw=lw, color=my_color)#, label=promoter_label[i])

# The spacer model
sb = bm.SpacerBinding(**{'k_on': 1.0, 'superhelical_op': -0.06, 'spacer': 17})
y = sb.rate_modulation(superhelical=x)
ax.plot(x, y, lw=lw, color=spacer_color, label='Spacer Model')

ax.text(-0.1, 1.1, olabel, transform=ax.transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

ax.grid(True)
ax.legend(loc='best',fontsize=font_size)
ax.set_xlabel('Superhelical Density', fontsize=xlabel_size)
ax.set_ylabel(r'Rate Modulation', fontsize=xlabel_size)
#ax.set_title('Binding Rate (Gaussian Modulation)', fontsize=title_size)
ax.set_title(mytitle, fontsize=title_size)
# c) Gaussian V2
# ----------------------------------------------------
V_num = 2
ax = axs[2,0]
olabel = outside_label[2]
mytitle = 'Gaussian Modulation: V2'
case = info_list[V_num]
responses_list_V2 = []
# Go through promoter
for j, pm in enumerate(promoter_cases):


    my_color = colors[j]

    figtitle = case['label'] + ': ' + promoter_labels[j] + ' Promoter'
    ax.set_title(figtitle, fontsize=title_size)  # color=colors[i],

    # Experimental susceptibility
    exp = pd.read_csv(susceptibility_files[0])  # read
    exp_x = exp['distance']
    exp_y = exp['Signal']
    exp_ys = exp['Error']

    # Load and prepare data, sort losses and find best case according to the weighted loss
    # --------------------------------------------------------------------
    with open(case[pm], 'rb') as file:
        data = pickle.load(file).results

    with open(case[pm], 'rb') as file:
        vals_dictlist = pickle.load(file).vals

    # Change std for error in the mean (or standard error)
    for dat in data:
        dat['prod_rate'][:, 1] = dat['prod_rate'][:, 1] / np.sqrt(60)

    # Filter data to include only entries with 'status' == 'ok'
    results = [item for item in data if item['status'] == 'ok']

    # Calculate new loss and just added to the list with results "results"
    for dat in results:
        new_loss = get_loss([dat], exp_y, exp_ys, w_mean=mean_weight, w_err=error_weight)
        dat['new_loss'] = new_loss[0]

    # Filter data_list to include only entries with 'status' == 'ok'
    dat_list = [item for item in results if item['status'] == 'ok']

    # Find the dictionary with the smallest 'loss'
    dat = min(dat_list, key=lambda x: x['new_loss'])

    # Get the index of this min_loss_data in the original data_list and append it
    min_loss_index = next(index for index, item in enumerate(dat_list) if item == dat)

    # Collect responses dict
    # ------------------------------------------------------------------------------
    # Let's form the response dict
    resp_dict = {}
    # Getting values from trials
    for key in vals_dictlist:
        resp_dict[key] = vals_dictlist[key][min_loss_index]

    # Getting values from promoter response
    for key in ['width', 'threshold']:
        resp_dict[key] = responses[j][key]
    # Rename kon to k_on
    resp_dict['k_on'] = resp_dict.pop('kon')

    # Let's save the responses
    responses_list_V2.append(resp_dict)

    # Plot Gaussian
    # ------------------------------------------------------------------------------
    # Gaussian binding
    gy = bm.GaussianBinding(**resp_dict).rate_modulation(superhelical=x)
    y =  (gy - np.min(gy)) / (np.max(gy) - np.min(gy)) # Let's normalize
    ax.plot(x, y, lw=lw, color=my_color)#, label=promoter_label[i])

# The spacer model
sb = bm.SpacerBinding(**{'k_on': 1.0, 'superhelical_op': -0.06, 'spacer': 17})
y = sb.rate_modulation(superhelical=x)
ax.plot(x, y, lw=lw, color=spacer_color, label='Spacer Model')

ax.text(-0.1, 1.1, olabel, transform=ax.transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

ax.grid(True)
#ax.legend(loc='best',fontsize=font_size)
ax.set_xlabel('Superhelical Density', fontsize=xlabel_size)
ax.set_ylabel(r'Rate Modulation', fontsize=xlabel_size)
ax.set_title(mytitle, fontsize=title_size)

# d) Bars V0
# ----------------------------------------------------
V_num = 0
ax = axs[0,1]
olabel = outside_label[3]
mytitle = 'Transition Rates: V0'
case = info_list[V_num]
responses_list_V0 = []
# Go through promoter and get respones
for j, pm in enumerate(promoter_cases):

    # Load and prepare data, sort losses and find best case according to the weighted loss
    # --------------------------------------------------------------------
    with open(case[pm], 'rb') as file:
        data = pickle.load(file).results

    with open(case[pm], 'rb') as file:
        vals_dictlist = pickle.load(file).vals

    # Change std for error in the mean (or standard error)
    for dat in data:
        dat['prod_rate'][:, 1] = dat['prod_rate'][:, 1] / np.sqrt(60)

    # Filter data to include only entries with 'status' == 'ok'
    results = [item for item in data if item['status'] == 'ok']

    # Calculate new loss and just added to the list with results "results"
    for dat in results:
        new_loss = get_loss([dat], exp_y, exp_ys, w_mean=mean_weight, w_err=error_weight)
        dat['new_loss'] = new_loss[0]

    # Filter data_list to include only entries with 'status' == 'ok'
    dat_list = [item for item in results if item['status'] == 'ok']

    # Find the dictionary with the smallest 'loss'
    dat = min(dat_list, key=lambda x: x['new_loss'])

    # Get the index of this min_loss_data in the original data_list and append it
    min_loss_index = next(index for index, item in enumerate(dat_list) if item == dat)

    # Collect responses dict
    # ------------------------------------------------------------------------------
    # Let's form the response dict
    resp_dict = {}
    # Getting values from trials
    for key in vals_dictlist:
        resp_dict[key] = vals_dictlist[key][min_loss_index]

    # Getting values from promoter response
    for key in ['width', 'threshold']:
        resp_dict[key] = responses[j][key]
    # Rename kon to k_on
    resp_dict['k_on'] = resp_dict.pop('kon')

    # Let's save the responses
    responses_list_V0.append(resp_dict)

# Sort the responses and plot
# -------------------------------
responses_list = responses_list_V0
# Remove the keys from all dictionaries in the list
for case_dict in responses_list:
    for key in nokeys:
        if key in case_dict:
            del case_dict[key]

keys = list(responses_list[0].keys())

# Number of cases and keys
n_cases = len(responses)
n_keys = len(keys)

# Set up the positions for each bar group (x-axis)
x = np.arange(n_keys)  # Position of each group on the x-axis
width = 0.8 / n_cases  # Dynamically calculate the width of each bar

# Plot each case in the list
for i, case_dict in enumerate(responses_list):

    values = [case_dict[key] for key in keys]  # Get values for the current case
    ax.bar(x + i * width - width * (n_cases - 1) / 2, values, width, label=promoter_cases[i], color=colors[i])

# Add labels, title, and custom ticks
ax.text(-0.1, 1.1, olabel, transform=ax.transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')
ax.set_xlabel('Rate Type', fontsize=xlabel_size)
ax.set_ylabel(r'Rate ($s^{-1}$)', fontsize=xlabel_size)
ax.set_title(mytitle, fontsize=title_size)
ax.set_xticks(x)
ax.set_xticklabels(keys)
#ax.set_ylim(t_rates_xlim)

# e) Bars V1 - 4
# ----------------------------------------------------
ax = axs[1,1]
olabel = outside_label[4]
responses_list = responses_list_V1
mytitle = 'Transition Rates: V1'

# Remove the keys from all dictionaries in the list
for case_dict in responses_list:
    for key in nokeys:
        if key in case_dict:
            del case_dict[key]

keys = list(responses_list[0].keys())

# Number of cases and keys
n_cases = len(responses)
n_keys = len(keys)

# Set up the positions for each bar group (x-axis)
x = np.arange(n_keys)  # Position of each group on the x-axis
width = 0.8 / n_cases  # Dynamically calculate the width of each bar

# Plot each case in the list
for i, case_dict in enumerate(responses_list):
    values = [case_dict[key] for key in keys]  # Get values for the current case
    ax.bar(x + i * width - width * (n_cases - 1) / 2, values, width, label=promoter_cases[i], color=colors[i])

# Add labels, title, and custom ticks
ax.text(-0.1, 1.1, olabel, transform=ax.transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')
ax.set_xlabel('Rate Type', fontsize=xlabel_size)
ax.set_ylabel(r'Rate ($s^{-1}$)', fontsize=xlabel_size)
ax.set_title(mytitle, fontsize=title_size)
ax.set_xticks(x)
ax.set_xticklabels(keys)
#ax.set_ylim(t_rates_xlim)

# f) Bars V2
# ----------------------------------------------------
ax = axs[2,1]
olabel = outside_label[5]
mytitle = 'Transition Rates: V2'

responses_list = responses_list_V2

# Remove the keys from all dictionaries in the list
for case_dict in responses_list:
    for key in nokeys:
        if key in case_dict:
            del case_dict[key]

keys = list(responses_list[0].keys())

# Number of cases and keys
n_cases = len(responses)
n_keys = len(keys)

# Set up the positions for each bar group (x-axis)
x = np.arange(n_keys)  # Position of each group on the x-axis
width = 0.8 / n_cases  # Dynamically calculate the width of each bar

# Plot each case in the list
for i, case_dict in enumerate(responses_list):
    values = [case_dict[key] for key in keys]  # Get values for the current case
    ax.bar(x + i * width - width * (n_cases - 1) / 2, values, width, label=promoter_cases[i], color=colors[i])


# Add labels, title, and custom ticks
ax.text(-0.1, 1.1, olabel, transform=ax.transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')
ax.set_xlabel('Rate Type', fontsize=xlabel_size)
ax.set_ylabel(r'Rate ($s^{-1}$)', fontsize=xlabel_size)
ax.set_title(mytitle, fontsize=title_size)
ax.set_xticks(x)
ax.set_xticklabels(keys)
#ax.set_ylim(t_rates_xlim)

plt.savefig(outfig+'.png')

plt.show()
