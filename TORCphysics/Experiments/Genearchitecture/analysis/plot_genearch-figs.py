import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns
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
outfile = 'define_new_loss'
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
for pcase in promoter_cases:
    susceptibility_files.append('../junier_data/'+ pcase + '.csv')


# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 4#5
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


# This one calculates one component per experiment. so basically, it is the mean distances and err distances (x,y)
def get_xy_loss_components(data_trials, exp_mean, exp_err, w_mean=1, w_err=1):

    # x and y are just the squared components, while xm and ym are the weighted ones.
    x_list = []
    y_list = []
    xm_list = []
    ym_list = []
    n = len(exp_mean)
    # First calculate susceptibility
    for trial in data_trials:
        pr = trial['prod_rate'][:,0]
        prs = trial['prod_rate'][:,1]

        y = pr / pr[4]  # susceptibility
        ys = np.zeros_like(y)

        for j in range(len(y)):
            ys[j] = y[j] * np.sqrt(np.square(prs[j] / pr[j]) + np.square(prs[4] / pr[4]))

        xloss = 0
        yloss = 0
        for i in range(n):
            xloss += np.square(exp_mean[i] - y[i])
            yloss += np.square(exp_err[i] - ys[i])

        x_list.append(xloss)
        y_list.append(yloss)
        xm_list.append(xloss*w_mean)
        ym_list.append(yloss*w_err)

    return np.array(x_list), np.array(y_list), np.array(xm_list), np.array(ym_list)

# This one does it for every single measurement, so you get like 12 datapoints per experiment
def get_xy_loss_components_all(data_trials, exp_mean, exp_err, w_mean=1, w_err=1):

    # x and y are just the squared components, while xm and ym are the weighted ones.
    x_list = []  # mean
    y_list = []  # err
    xm_list = []
    ym_list = []
    n = len(exp_mean)
    # First calculate susceptibility
    for trial in data_trials:
        pr = trial['prod_rate'][:,0]
        prs = trial['prod_rate'][:,1]

        y = pr / pr[4]  # susceptibility
        ys = np.zeros_like(y)

        for j in range(len(y)):
            ys[j] = y[j] * np.sqrt(np.square(prs[j] / pr[j]) + np.square(prs[4] / pr[4]))

        # Now, let's calculate the x and y components
        for i in range(n):
            a = np.square(exp_mean[i] - y[i])
            b = np.square(exp_err[i] - ys[i])
            x_list.append(a)
            y_list.append(b)
            xm_list.append(a*w_mean)
            ym_list.append(b*w_err)
    return np.array(x_list), np.array(y_list), np.array(xm_list), np.array(ym_list)

# Gets a quantity x for all values
def get_xquantity_all(data_trials, xquantity):
    x_list = []
    for trial in data_trials:
        n = len(trial[xquantity])
        for i in range(n):
            x_list.append(trial[xquantity][i])
    return np.array(x_list)


# Process and plot - Let's plot susceptibility vs upstream barrier distance of best cases
#-----------------------------------------------------------------------------------------------------------------------
outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)']
ncols = len(info_list)
nrows = len(promoter_cases) # 1.- loss vs prodrate, 2.- xy weighted components, 3.- xy not weighted,
# 4.- prod rate vs susceptibility, 5.- prod rate vs susceptibility mean, 6.- prod rate vs susceptibilty, not mean
outfig = 'genearch_susceptibility'
fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*width, nrows*height),
tight_layout=True, sharex=True, sharey=True)

# Go through case - V0, V1, V2
# ---------------------------------------------------
s= -1
for i, case in enumerate(info_list):

    # Go through promoter
    for j, pm in enumerate(promoter_cases):

        s=s+1 # This is for the outside label

        ax = axs[j, i]
        my_color = colors[j]

        figtitle = case['label'] + ': ' + promoter_labels[j] + ' Promoter'
        ax.set_title(figtitle, fontsize=title_size)  # color=colors[i],

        # Experimental susceptibility
        exp = pd.read_csv(susceptibility_files[j])  # read
        exp_x = exp['distance']
        exp_y = exp['Signal']
        exp_ys = exp['Error']


        # Load and prepare data
        # --------------------------------------------------------------------
        with open(case[pm], 'rb') as file:
            data = pickle.load(file).results

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

        # Plot
        x = dat['distance']
        pr = dat['prod_rate'][:, 0]
        prs = dat['prod_rate'][:, 1]  # /np.sqrt(60)
        y = pr / pr[4]  # susceptibility
        ys = np.zeros_like(y)
        for j in range(len(y)):
            ys[j] = y[j] * np.sqrt(np.square(prs[j] / pr[j]) + np.square(prs[4] / pr[4]))

        ax.plot(x, y, model_ls, lw=lw, color=my_color, label='TORCPhysics')
        ax.fill_between(x, y - ys, y + ys, lw=lw, color=my_color, alpha=0.2)

        # Load experimental and plot
        x = exp_x
        y = exp_y
        ys = exp_ys
        ax.errorbar(x, y, yerr=ys, fmt=exp_ls, capsize=5, color='black', label='Experiment')

        # Loss
        error = dat['new_loss']
        formatted_sum_err = "{:.3g}".format(error)  # Formatting to three significant figures
        er_label = r'$\epsilon={}$'.format(formatted_sum_err)

        # Add the text box to the plot
        props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
        ax.text(0.7, 0.1, er_label, transform=ax.transAxes, fontsize=font_size,
                verticalalignment='top', bbox=props)

        # Add label outside the plot
        ax.text(-0.05, 1.1, outside_label[s], transform=ax.transAxes,
                fontsize=font_size * 1.5, fontweight='bold', va='center', ha='center')

        ax.grid(True)
        ax.set_xlim([80, 6000])
        ax.set_xscale('log')
for i in range(3):
    axs[i,0].set_ylabel(r'Susceptibility', fontsize=xlabel_size)
    axs[2,i].set_xlabel('Upstream Distance (bp)', fontsize=xlabel_size)
axs[0,0].legend(loc='best')

plt.savefig(outfig+'.png')
plt.show()

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Process and plot - Prod rate vs loss
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)']
ncols = len(info_list)
nrows = len(promoter_cases) # 1.- loss vs prodrate, 2.- xy weighted components, 3.- xy not weighted,
# 4.- prod rate vs susceptibility, 5.- prod rate vs susceptibility mean, 6.- prod rate vs susceptibilty, not mean
outfig = 'genearch_prod-rate_vs_loss'
fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*width, nrows*height),
tight_layout=True, sharex=True, sharey=True)


# Go through case - V0, V1, V2
# ---------------------------------------------------
s= -1
for i, case in enumerate(info_list):

    # Go through promoter
    for j, pm in enumerate(promoter_cases):

        s=s+1 # This is for the outside label

        ax = axs[j, i]
        my_color = colors[j]

        figtitle = case['label'] + ': ' + promoter_labels[j] + ' Promoter'
        ax.set_title(figtitle, fontsize=title_size)  # color=colors[i],

        # Load and prepare data
        # --------------------------------------------------------------------
        with open(case[pm], 'rb') as file:
            data = pickle.load(file).results

        # Change std for error in the mean (or standard error)
        for dat in data:
            dat['prod_rate'][:, 1] = dat['prod_rate'][:, 1] / np.sqrt(60)

        # Filter data to include only entries with 'status' == 'ok'
        results = [item for item in data if item['status'] == 'ok']

        # Calculate new loss and just added to the list with results "results"
        for dat in results:
            new_loss = get_loss([dat], exp_y, exp_ys, w_mean=mean_weight, w_err=error_weight)
            dat['new_loss'] = new_loss[0]

        # ---------------------------------------------------------------------------------
        # Collect variables
        # ---------------------------------------------------------------------------------
        # Adjust number of considered cases so all dists have the same number of values
        n = len(data)
        nconsidered = int(n * percentage_threshold)

        # Collect variables
        loss = get_loss(results, exp_y, exp_ys, w_mean=mean_weight, w_err=error_weight)
        prod_rate = [np.mean(t['prod_rate'][:, 0]) for t in results]
        local_superhelical = [np.mean(t['local_superhelical']) for t in results]
        global_superhelical = [np.mean(t['global_superhelical']) for t in results]
        nRNAPs = [np.mean(t['nenzymes']['RNAP']) for t in results]
        ntopoIs = [np.mean(t['nenzymes']['topoI']) for t in results]
        ngyrases = [np.mean(t['nenzymes']['gyrase']) for t in results]
        xloss, yloss, xwloss, ywloss = get_xy_loss_components(results, exp_y, exp_ys, w_mean=mean_weight,
                                                              w_err=error_weight)
        susceptibility = [np.mean(t['susceptibility'][:, 0]) for t in results]

        # And put them in df
        case_df = pd.DataFrame({
            "loss": loss,
            "prod_rate": prod_rate,
        })

        # ---------------------------------------------------------------------------------
        # Filter according loss
        # ---------------------------------------------------------------------------------
        # First get rid of NaNs
        case_df = case_df.dropna(subset=['loss'])
        # Then filter
        case_df = case_df.sort_values(by='loss', ascending=False)  # , inplace=True)
        err_threshold = case_df['loss'].iloc[-nconsidered]
        # Filter according error
        filtered_df = case_df[case_df['loss'] <= err_threshold]

        # ----------------------------------------------------
        ylabel = 'prod_rate'
        xlabel = 'loss'

        x = case_df[xlabel].to_numpy()
        y = case_df[ylabel].to_numpy()

        xf = filtered_df[xlabel].to_numpy()
        yf = filtered_df[ylabel].to_numpy()

        ax.plot(x, y, 'o', color=my_color, alpha=0.3, label='cases')
        #ax.plot(x, y, 'o', color='gray', alpha=0.3, label='cases')
        #ax.plot(xf, yf, 'o', color=my_color, alpha=0.3, label='filtered')
        ax.set_xscale('log')
        ax.grid(True)

        # Add label outside the plot
        ax.text(-0.1, 1.15, outside_label[s], transform=ax.transAxes,
                fontsize=font_size * 1.5, fontweight='bold', va='center', ha='center')

for i in range(3):
    axs[i,0].set_ylabel(r'Production Rate ($s^{-1}$)', fontsize=xlabel_size)
    axs[2,i].set_xlabel('Loss', fontsize=xlabel_size)
axs[0,0].legend(loc='best')

plt.savefig(outfig+'.png')
plt.show()

# Process and plot - Let's plot number of enzymes as a function of upstream barrier
#-----------------------------------------------------------------------------------------------------------------------
outfig = 'genearch_nenzymes'
fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*width, nrows*height),
tight_layout=True, sharex=True, sharey=True)

enzymes_names = ['gyrase', 'topoI', 'RNAP']

enzyme_colors = ['blue', 'red', 'black']

# Go through case - V0, V1, V2
# ---------------------------------------------------
s= -1
for i, case in enumerate(info_list):

    # Go through promoter
    for j, pm in enumerate(promoter_cases):

        s=s+1 # This is for the outside label

        ax = axs[j, i]

        figtitle = case['label'] + ': ' + promoter_labels[j] + ' Promoter'
        ax.set_title(figtitle, fontsize=title_size)  # color=colors[i],

        # Load and prepare data
        # --------------------------------------------------------------------
        with open(case[pm], 'rb') as file:
            data = pickle.load(file).results

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

        # Plot
        x = dat['distance']
        nenzymes = dat['nenzymes']

        for k, ename in enumerate(enzymes_names):
            y = nenzymes[ename][:,0]
            ys = nenzymes[ename][:,1]
            my_color = enzyme_colors[k]

            # ax.errorbar(x, y, yerr=ys, fmt=exp_ls, capsize=5, color=my_color, label=ename)

            ax.plot(x, y, model_ls, lw=lw, color=my_color, label=ename)
            if ename=='RNAP':
                ax.fill_between(x, y - ys, y + ys, lw=lw, color=my_color, alpha=0.2)

        # Add label outside the plot
        ax.text(-0.05, 1.15, outside_label[s], transform=ax.transAxes,
                fontsize=font_size * 1.5, fontweight='bold', va='center', ha='center')

        ax.grid(True)
        ax.set_xlim([80, 6000])
        ax.set_xscale('log')
for i in range(3):
    axs[i,0].set_ylabel(r'Bound Enzymes', fontsize=xlabel_size)
    axs[2,i].set_xlabel('Upstream Distance (bp)', fontsize=xlabel_size)
axs[0,0].legend(loc='best')

plt.savefig(outfig+'.png')
plt.show()


# Process and plot - Let's plot number of enzymes as a function of upstream barrier
#-----------------------------------------------------------------------------------------------------------------------
outfig = 'genearch_superhelical'
fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*width, nrows*height),
tight_layout=True, sharex=True, sharey=True)

# Go through case - V0, V1, V2
# ---------------------------------------------------
s= -1
for i, case in enumerate(info_list):

    # Go through promoter
    for j, pm in enumerate(promoter_cases):

        s=s+1 # This is for the outside label

        ax = axs[j, i]

        figtitle = case['label'] + ': ' + promoter_labels[j] + ' Promoter'
        ax.set_title(figtitle, fontsize=title_size)  # color=colors[i],

        # Load and prepare data
        # --------------------------------------------------------------------
        with open(case[pm], 'rb') as file:
            data = pickle.load(file).results

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

        # Plot
        x = dat['distance']

        y = dat['global_superhelical'][:,0]
        ys = dat['global_superhelical'][:,1]
        ax.errorbar(x, y, yerr=ys, fmt=model_ls, capsize=5, color='black', label='global')

        y = dat['local_superhelical'][:,0]
        ys = dat['local_superhelical'][:,1]
        ax.errorbar(x, y, yerr=ys, fmt=model_ls, capsize=5, color='red', label='local')

        #ax.plot(x, y, model_ls, lw=lw, color=my_color, label=ename)
        #ax.fill_between(x, y - ys, y + ys, lw=lw, color=my_color, alpha=0.2)

        # Add label outside the plot
        ax.text(-0.05, 1.15, outside_label[s], transform=ax.transAxes,
                fontsize=font_size * 1.5, fontweight='bold', va='center', ha='center')

        ax.grid(True)
        ax.set_xlim([80, 6000])
        ax.set_ylim([-.1, 0])
        ax.set_xscale('log')
for i in range(3):
    axs[i,0].set_ylabel(r'Superhelical Density', fontsize=xlabel_size)
    axs[2,i].set_xlabel('Upstream Distance (bp)', fontsize=xlabel_size)
axs[0,0].legend(loc='best')

plt.savefig(outfig+'.png')
plt.show()
