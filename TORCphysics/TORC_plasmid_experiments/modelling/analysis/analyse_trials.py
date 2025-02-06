import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns

# Description
#-----------------------------------------------------------------------------------------------------------------------
# I want to plot bar plots of the parameters inferred

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
dir_source = '../optimization/third_strategy.2/'
# Nuevos
v_code = 'block-full-dist_op_TORC_plasmid_st3.2-02'
v_code = 'block-full-trackingON-dist_op_TORC_plasmid_st3.2-02'

# Viejos
#v_code = 'block-dist_op_TORC_plasmid_v2-01'
#v_code = 'block-dist_op_TORC_plasmid_v2-02'
#v_code = 'full-dist_op_TORC_plasmid_st3-02'
#v_code = 'full-trackingON-dist_op_TORC_plasmid_st3-02'

trials_file = dir_source + v_code + '-trials.pkl'
params_file = dir_source + v_code + '.csv'
loss_file = dir_source + v_code + '-values.csv'
#ref_pkl_file  = dir_source + 'full-dist_op_TORC_plasmid_st3-01'+ '.pkl' # We just want a pkl file to load the reference dist
pkl_file  = dir_source + v_code+ '.pkl' # We just want a pkl file to load the reference dist

out_file = 'trials_'+v_code

percentage_threshold = .05
#err_threshold = 0.5 # The minimum error we want


# Load
#-----------------------------------------------------------------------------------------------------------------------
with open(trials_file, 'rb') as file:
    trials_data = pickle.load(file)

with open(pkl_file, 'rb') as file:
    pkl_data = pickle.load(file)

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 8
height = 5
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16


colors = ['green', 'blue', 'red']


#-----------------------------------------------------------------------------------------------------------------------
# PROCESS
#-----------------------------------------------------------------------------------------------------------------------

# Let's sort the losses
results = trials_data.results

system_loss_df = pd.DataFrame([t['system_loss'] for t in results])

loss_df = pd.DataFrame({'loss':[t['loss'] for t in results]})

# Assuming loss_df is a single-column DataFrame
system_loss_df['loss'] = loss_df.squeeze()  # Convert loss_df to Series if needed

system_loss_df = system_loss_df.sort_values(by='loss', ascending=False)#, inplace=True)

n = len(system_loss_df['loss'])
nconsidered = int(n*percentage_threshold)
err_threshold = system_loss_df['loss'].iloc[-nconsidered]
print('Number of tests', n)
print('Considered', nconsidered)
print('For ', percentage_threshold*100, '%')
# Filter according error
filtered_df = system_loss_df[system_loss_df['loss'] <= err_threshold]

# Set with minimum loss
dat = min(results, key=lambda x: x['loss'])

# Find the index of the dictionary with minimum loss
index = next(i for i, item in enumerate(results) if item == dat)


# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we load
fig, axs = plt.subplots(5, figsize=(width, 5*height), tight_layout=True)

# Loss
# ----------------------------------------------------------------------------------------------------------------------
ms=6
ax = axs[0]
ax.set_title('loss for '+ v_code)

x = np.arange(1, n+1, 1)
loss = system_loss_df['loss'].to_numpy()
ax.plot(x,loss, 'o', ms=ms, color='blue')
ax.plot(x[n-nconsidered:], loss[n-nconsidered:], 'o', ms=ms, color='red')

ax.grid(True)
ax.set_xlabel('test')
ax.set_ylabel('loss')
#ax.set_yscale('log')

# System loss distribution
# ----------------------------------------------------------------------------------------------------------------------
ax =axs[1]
ax.set_title('System loss for '+ v_code)

system_loss = system_loss_df.drop('loss', axis=1)

sns.violinplot(data=system_loss, ax=ax, inner="quart")#, cut=0, color=colors[i])
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_ylabel('System loss')
ax.grid(True)

# Filtered system loss distribution
# ----------------------------------------------------------------------------------------------------------------------
ax =axs[2]
ax.set_title('Filtered loss for '+ v_code)

system_loss = filtered_df.drop('loss', axis=1)

sns.violinplot(data=system_loss, ax=ax, inner="quart")#, cut=0, color=colors[i])
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_ylabel('System loss')
ax.grid(True)

# Relatives - Comparing with experiments
# ----------------------------------------------------------------------------------------------------------------------
ax = axs[3]
#keys = list(results[0]['system_loss'].keys())
#relative_rates = trials_data.results[-1]['relative_rate']
keys = [my_dict['name'] for my_dict in pkl_data]
print('keys', keys)

# Number of cases and keys
n_cases = 2#len(pkl_data) # because we have 2 measurements per case - experimental and simulation
n_keys = len(keys)

# Set up the positions for each bar group (x-axis)
x = np.arange(n_keys)  # Position of each group on the x-axis
width = .8 / n_cases  # Dynamically calculate the width of each bar

labels = ['exp', 'sim']
colors=['red', 'blue']
# Plot each case in the list
for i, measurement in enumerate(['reference', 'relative_rate']):

    values = [my_dict[measurement] for my_dict in pkl_data] # mean
    std = np.std(values, axis=1)
    values = np.mean(values, axis=1)

    ax.bar(x + i * width - width * (n_cases - 1) / 2, values, width, yerr=std, label=labels[i], color=colors[i])

#    if measurement=='reference':
#        values = []
#        for name in keys:
##            for my_dict in ref_pkl_data:
#                if my_dict['name'] == name:
#                    values.append(my_dict['reference'])
##        std = np.std(values, axis=1)
#        values = np.mean(values, axis=1)
#    elif measurement=='relative_rate':
#        values = []
#        for name in keys:
#            values.append(relative_rates[name])
#        std = np.std(values, axis=1)
#        values = np.mean(relative_rates, axis=1)


#    ax.bar(x + i * width - width * (n_cases - 1) / 2, values, width, yerr=std, label=labels[i], color=colors[i])

# Let's draw the error
#objective = sum([my_dict['objective'] for my_dict in pkl_data])
# Add label outside the plot
# Add errors
#formatted_sum_err = "{:.3f}".format(objective)  # Formatting to three significant figures
#er_label = r'$\epsilon={}$'.format(formatted_sum_err)
#props = dict(boxstyle='round', facecolor='silver', alpha=0.4)
#axs.text(0.3, 0.95, er_label, transform=axs.transAxes, fontsize=font_size,
#        verticalalignment='top', bbox=props)
#print(objective)

# Add labels, title, and custom ticks
ax.set_xlabel('System')
ax.set_ylabel('Relative expression')
ax.set_title('TORC Plasmid')
ax.set_xticks(x)
ax.set_xticklabels(keys, rotation=45, ha='right')
ax.legend(loc='best')
ax.grid(True)
ax.set_ylim(0,2)

# Production rates - Comparing with experiments
# ----------------------------------------------------------------------------------------------------------------------
ax = axs[4]
keys = [my_dict['name'] for my_dict in pkl_data]
colors=['yellow', 'blue', 'red', 'green', 'purple']
gene_names = ['PleuWT', 'tetA', 'antitet', 'bla']#, 'lacI']

# Number of cases and keys
n_cases = len(gene_names) # because we have 2 measurements per case - experimental and simulation
n_keys = len(keys)

# Set up the positions for each bar group (x-axis)
x = np.arange(n_keys)  # Position of each group on the x-axis
bwidth = .8 / n_cases  # Dynamically calculate the width of each bar

labels = gene_names
# Plot each case in the list
for i, name in enumerate(gene_names):

    mean = [my_dict['production_rate'][name][0] for my_dict in pkl_data] # mean
    std = [my_dict['production_rate'][name][1] for my_dict in pkl_data] # std

    print(name)
    ax.bar(x + i * bwidth - bwidth * (n_cases - 1) / 2, mean, bwidth, yerr=std, label=labels[i], color=colors[i], alpha=0.7)


# Add labels, title, and custom ticks
ax.set_xlabel('System')
ax.set_ylabel(r'Production rate ($s^{-1}$)')
ax.set_title('TORC Plasmid')
ax.set_xticks(x)
ax.set_xticklabels(keys, rotation=45, ha='right')
ax.legend(loc='best')
ax.grid(True)

#plt.savefig(out_file+'.png')
#plt.savefig(out_file+'.pdf')
plt.show()
