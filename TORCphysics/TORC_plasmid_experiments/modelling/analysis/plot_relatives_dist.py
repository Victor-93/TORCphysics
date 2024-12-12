import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle


# Description
#-----------------------------------------------------------------------------------------------------------------------
# We want to plot the results from the TORC plasmid following the optmization process

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
dir_source = '../optimization/'
#v_code = 'block-dist_op_TORC_plasmid'
#v_code = 'batch-dist_op_TORC_plasmid-01'
#v_code = 'block-batch-dist_op_TORC_plasmid-01'
v_code = 'block-batch-dist_op_TORC_plasmid-04'
pkl_file  = dir_source + v_code+ '.pkl'
params_file = dir_source + v_code + '.csv'
loss_file = dir_source + v_code + '-values.csv'

#pkl_file  = dir_source + 'rep-TORC_plasmid.pkl'
#params_file = dir_source + 'calibration_TORC_plasmid.csv'
#loss_file = dir_source + 'calibration_TORC_plasmid-values.csv'

out_file = 'relatives_'+v_code

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 8
height = 5
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16


colors = ['green', 'blue', 'red']

# Load
#-----------------------------------------------------------------------------------------------------------------------
with open(pkl_file, 'rb') as file:
    pkl_data = pickle.load(file)


# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we do the process
fig, axs = plt.subplots(1, figsize=(width, 1*height), tight_layout=True, sharex=True)

keys = [my_dict['name'] for my_dict in pkl_data]

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

    axs.bar(x + i * width - width * (n_cases - 1) / 2, values, width, yerr=std, label=labels[i], color=colors[i])

# Let's draw the error
objective = sum([my_dict['objective'] for my_dict in pkl_data])
# Add label outside the plot
# Add errors
formatted_sum_err = "{:.3f}".format(objective)  # Formatting to three significant figures
er_label = r'$\epsilon={}$'.format(formatted_sum_err)
props = dict(boxstyle='round', facecolor='silver', alpha=0.4)
#axs.text(0.3, 0.95, er_label, transform=axs.transAxes, fontsize=font_size,
#        verticalalignment='top', bbox=props)
print(objective)

# Add labels, title, and custom ticks
axs.set_xlabel('System')
axs.set_ylabel('Relative expression')
axs.set_title('TORC Plasmid')
axs.set_xticks(x)
axs.set_xticklabels(keys, rotation=45, ha='right')
axs.legend(loc='best')
axs.grid(True)

#plt.savefig(out_file+'.png')
plt.show()



