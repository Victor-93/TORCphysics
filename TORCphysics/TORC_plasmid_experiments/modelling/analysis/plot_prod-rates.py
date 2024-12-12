import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle


# Description
#-----------------------------------------------------------------------------------------------------------------------
# We want to plot the expression rates from the TORC plasmid following the optmization process

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
dir_source = '../optimization/'
v_code = 'block-op_TORC_plasmid'
#v_code = 'block-rep-TORC_plasmid'
#v_code = 'op_TORC_plasmid'
v_code = 'block-dist_op_TORC_plasmid'
v_code = 'batch-dist_op_TORC_plasmid-01'


pkl_file  = dir_source + v_code+ '.pkl'

out_file = 'prod-rate_'+v_code

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 12#8
height = 5
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

colors=['yellow', 'blue', 'red', 'green', 'purple']
gene_names = ['PleuWT', 'tetA', 'antitet', 'bla']#, 'lacI']

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
    axs.bar(x + i * bwidth - bwidth * (n_cases - 1) / 2, mean, bwidth, yerr=std, label=labels[i], color=colors[i], alpha=0.7)


# Add labels, title, and custom ticks
axs.set_xlabel('System')
axs.set_ylabel(r'Production rate ($s^{-1}$)')
axs.set_title('TORC Plasmid')
axs.set_xticks(x)
axs.set_xticklabels(keys, rotation=45, ha='right')
axs.legend(loc='best')
axs.grid(True)
#plt.savefig(out_file+'.png')
plt.show()



