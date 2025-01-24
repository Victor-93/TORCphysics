import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
from TORCphysics import utils
from TORCphysics import binding_model as bm

# Description
#-----------------------------------------------------------------------------------------------------------------------
# Plots general plot for the best case of the three promoters. It uses the trials pkl.
# TODO: Do a special version for V0

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
dt=1.0#0.5

# Selec this for V2
path = 'genearch_V2/'
model_code='sus_genearch_V2-01-'
figtitle = model_code
outfile = model_code+'trials'

# Selec this for V2 - small dist
#path = 'genearch_V2_small-distance/'
#model_code='sus_genearch_V2-01-small-dist-'
#figtitle = model_code
#outfile = model_code+'trials'

# Select this one for V1
path = 'genearch_V1/'
model_code='sus_genearch_V1-01-'
figtitle = model_code + ' gamma0.1'
outfile = model_code+'trials'

# Selec this for V0
#path = 'genearch_V0/'
#model_code='sus_genearch_V0-01-'
#figtitle = model_code + ' gamma0.1'
#outfile = model_code+'trials'

susceptibility_files = []
trials_file = []
promoter_responses_files = []
for pcase in promoter_cases:
    susceptibility_files.append('../junier_data/'+ pcase + '.csv')

    # This is for V2
    #trials_file.append(path+model_code+pcase+'_dt'+str(dt)+'-trials.pkl')

    # This for V1 with gamma
    trials_file.append(path+model_code+pcase+'_gamma0.1_dt'+str(dt)+'-trials.pkl')

    # This for V0 with gamma
    #trials_file.append(path+model_code+pcase+'_gamma0.1_dt'+str(dt)+'-trials.pkl')

    # This contains the promoter responses
    promoter_responses_files.append('../promoter_responses/' + pcase + '.csv')


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
def open_rate(superhelical, k_open, threshold, width):
    U = utils.opening_function(superhelical, threshold, width)
    rate = k_open * np.exp(-U)
    return rate

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

# Load
#-----------------------------------------------------------------------------------------------------------------------
pickle_data = []
for pickle_file in trials_file:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_data.append(data)

responses = []
for file in promoter_responses_files:
    data = pd.read_csv(file)
    data = data.to_dict(orient='records')[0]  # Get the first (and only) dictionary
    responses.append(data)

# Plot
#-----------------------------------------------------------------------------------------------------------------------
sigma = np.arange(-.13, .0, 0.001)
# Let's plot as we do the process
fig, axs = plt.subplots(3,2, figsize=(2*width, 3*height), tight_layout=True)
fig.suptitle(figtitle, fontsize=title_size)

outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)']

# Plot susceptibility
# ----------------------------------------------------------------------------------
n = len(trials_file)
min_loss_index = []
for i in range(n):

    # Load experimental and plot
    exp = pd.read_csv(susceptibility_files[i]) # read

    ax = axs[i,0]

    ax.set_title(titles[i], fontsize=title_size) #color=colors[i],

    # Load data
    # ----------------------------------------------------------------------------------
    data_list = pickle_data[i].results

    # Calculate new loss
    for data in data_list:
        data['prod_rate'][:, 1] = data['prod_rate'][:, 1]/np.sqrt(60)  #Change standard deviation for standard error
        new_loss = get_loss([data], exp['Signal'], exp['Error'], w_mean=1.0, w_err=1.0)
        data['new_loss'] = new_loss[0]

    # Filter data_list to include only entries with 'status' == 'ok'
    data_list = [item for item in data_list if item['status'] == 'ok']

    # Find the dictionary with the smallest 'loss'
    #data = min(data_list, key=lambda x: x['loss'])
    data = min(data_list, key=lambda x: x['new_loss'])

    # Get the index of this min_loss_data in the original data_list and append it
    min_loss_index.append( next(index for index, item in enumerate(data_list) if item == data) )

    # Plot
    x = data['distance']
    pr = data['prod_rate'][:, 0]
    prs = data['prod_rate'][:, 1]#/np.sqrt(60)
    y = pr / pr[4]  # susceptibility
    ys = np.zeros_like(y)
    for j in range(len(y)):
        ys[j] = y[j] * np.sqrt(np.square(prs[j] / pr[j]) + np.square(prs[4] / pr[4]))

    if i == 0:
        ax.plot(x, y, model_ls, lw=lw, color=colors[i], label='TORCPhysics')
    else:
        ax.plot(x, y, model_ls, lw=lw, color=colors[i])
    ax.fill_between(x, y-ys, y+ys, lw=lw, color=colors[i], alpha=0.2)


    # Load experimental and plot
    # exp = pd.read_csv(susceptibility_files[i]) # read
    x = exp['distance']
    y = exp['Signal']
    ys = exp['Error']
    if i == 0:
        ax.errorbar(x, y, yerr=ys, fmt=exp_ls, capsize=5,color=colors[i], label='Experiment')
    else:
        ax.errorbar(x, y, yerr=ys, fmt=exp_ls, capsize=5,color=colors[i])

    # Error quantification - meansquared error
    error = np.sum((exp['Signal'] - data['susceptibility'][:,0]) ** 2)/len(y)
    # Loss
    error = data['new_loss']
    formatted_sum_err = "{:.3g}".format(error)  # Formatting to three significant figures
    er_label = r'$\epsilon={}$'.format(formatted_sum_err)

    # Add the text box to the plot
    props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
    ax.text(0.7, 0.1, er_label, transform=ax.transAxes, fontsize=font_size,
            verticalalignment='top', bbox=props)

    # Add label outside the plot
    ax.text(-0.1, 1.1, outside_label[i], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

    ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    ax.set_ylabel(r'Susceptibility', fontsize=font_size)
    ax.grid(True)
    ax.set_xlim([80, 6000])
    ax.set_xscale('log')

axs[0,0].legend(loc='best')

# Responses part
#-----------------------------------------------------------------------------------------------------------------------
x = sigma
promoter_label = ['Weak', 'Medium', 'Strong']
# Let's first load the responses to a dict
responses_list = []
for i in range(n):

    # Let's form the response dict
    resp_dict = {}
    # Getting values from trials
    vals_dictlist = pickle_data[i].vals
    for key in vals_dictlist:
        resp_dict[key] = vals_dictlist[key][min_loss_index[i]]
    # Getting values from promoter response
    for key in ['width', 'threshold']:
        resp_dict[key] = responses[i][key]
    # Rename kon to k_on
    resp_dict['k_on'] = resp_dict.pop('kon')

    responses_list.append(resp_dict)

# Now let's plot the responses
for i in range(n):

    resp_dict = responses_list[i]
    # Gaussian binding
    gy = bm.GaussianBinding(**resp_dict).rate_modulation(superhelical=x)
    y =  (gy - np.min(gy)) / (np.max(gy) - np.min(gy)) # Let's normalize
    axs[0,1].plot(x, y, lw=lw, color=colors[i])#, label=promoter_label[i])

    # Openning rate
    oy = open_rate(x, resp_dict['k_open'], resp_dict['threshold'], resp_dict['width'])
    y =  (oy - np.min(oy)) / (np.max(oy) - np.min(oy)) # Let's normalize
    axs[1,1].plot(x, y, lw=lw, color=colors[i], label=promoter_label[i])

# The spacer model
sb = bm.SpacerBinding(**{'k_on':1.0, 'superhelical_op': -0.06, 'spacer': 17})
y = sb.rate_modulation(superhelical=x)
axs[0,1].plot(x, y, lw=lw, color='gray', label='Spacer Model')

# Meyer model for PelE? model
mb= bm.MeyerPromoterOpening (**{'k_on':1.0})
y = mb.rate_modulation(superhelical=x)
y = (y - np.min(y)) / (np.max(y) - np.min(y))  # Let's normalize
axs[1,1].plot(x, y, lw=lw, color='gray', label='pelE')

axs[0,1].text(-0.1, 1.1, outside_label[3], transform=axs[0,1].transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')
axs[1,1].text(-0.1, 1.1, outside_label[4], transform=axs[1,1].transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

axs[0,1].grid(True)
axs[1,1].grid(True)
axs[0,1].legend(loc='best',fontsize=font_size)
axs[1,1].legend(loc='best',fontsize=font_size)
axs[0,1].set_xlabel('Superhelical Density', fontsize=xlabel_size)
axs[1,1].set_xlabel('Superhelical Density', fontsize=xlabel_size)
axs[0,1].set_ylabel(r'Rate Modulation', fontsize=xlabel_size)
axs[1,1].set_ylabel(r'Rate Modulation', fontsize=xlabel_size)
axs[0,1].set_title('Binding Rate (Gaussian Modulation)', fontsize=title_size)
axs[1,1].set_title('Open Complex Formation Rate (SIST)', fontsize=title_size)

# Bars of values
#-----------------------------------------------------------------------------------------------------------------------
ax = axs[2,1]

nokeys = ['superhelical_op', 'spread', 'width', 'threshold', 'gamma', 'kappa', 'velocity', 'stall_torque']
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
ax.set_xlabel('Rate Type', fontsize=xlabel_size)
ax.set_ylabel(r'Rate ($s^{-1}$)', fontsize=xlabel_size)
ax.set_title('Transition Rates', fontsize=title_size)
ax.set_xticks(x)
ax.set_xticklabels(keys)

axs[2,1].text(-0.1, 1.1, outside_label[5], transform=axs[2,1].transAxes,
        fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

#plt.savefig(outfile+'.png')
plt.show()
