import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle


# Description
#-----------------------------------------------------------------------------------------------------------------------
# I want to plot bar plots of the parameters inferred

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
dir_source = '../optimization/'
v_code = 'block-op_TORC_plasmid'
pkl_file  = dir_source + v_code+ '.pkl'
params_file = dir_source + v_code + '.csv'
loss_file = dir_source + v_code + '-values.csv'

out_file = 'TORC_plasmid-test'

err_threshold = 0.3 # The minimum error we want

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

# Read inputs
#-----------------------------------------------------------------------------------------------------------------------
df = pd.read_csv(loss_file)
#df = df.sort_values(by='loss', ascending=False)#, inplace=True)
# Filter according error
filtered_df = df[df['loss'] <= err_threshold]

# Calculate averages and standard deviations
df_avg = pd.DataFrame(filtered_df.mean(axis=0))
df_std = pd.DataFrame(filtered_df.std(axis=0))

df_avg = filtered_df.mean(axis=0).to_frame().T
df_std = filtered_df.std(axis=0).to_frame().T

# Let's sort the rates
promoter_keys = ['k_on', 'k_off', 'k_open', 'k_closed',  'k_ini', 'bridge_on', 'bridge_off']
gene_names = ['PleuWT', 'tetA', 'antitet', 'bla']#, 'lacI']

# Add new keys to the genes, but they are 0 because they don't have that parameter
for name in gene_names:
    df_avg[name+'_bridge_on'] = 0.0
    df_std[name + '_bridge_on'] = 0.0
    df_avg[name+'_bridge_off'] = 0.0
    df_std[name + '_bridge_off'] = 0.0
gene_names.append('lacI') # Let's add lacI - even though is not a gene
# Now add k_closed, k_open, and 'k_ini to lacI
df_avg['lacI_k_open'] = 0.0
df_avg['lacI_k_closed'] = 0.0
df_avg['lacI_k_ini'] = 0.0
df_std['lacI_k_open'] = 0.0
df_std['lacI_k_closed'] = 0.0
df_std['lacI_k_ini'] = 0.0
df_avg.rename(columns={'bridge_on': 'lacI_bridge_on'}, inplace=True)
df_avg.rename(columns={'bridge_off': 'lacI_bridge_off'}, inplace=True)
df_std.rename(columns={'bridge_on': 'lacI_bridge_on'}, inplace=True)
df_std.rename(columns={'bridge_off': 'lacI_bridge_off'}, inplace=True)

# Number of cases and keys
n_cases = len(gene_names) # because we have 2 measurements per case - experimental and simulation
n_keys = len(promoter_keys)

# Set up the positions for each bar group (x-axis)
x = np.arange(n_keys)  # Position of each group on the x-axis
bwidth = .85 / n_cases  # Dynamically calculate the width of each bar

labels = gene_names
colors=['yellow', 'blue', 'red', 'green', 'purple']

# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we load
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True, sharex=True)

for i, name in enumerate(gene_names):

    mean = []
    for key in promoter_keys:
        mean.append( df_avg[name + '_' + key][0])
    std = []
    for key in promoter_keys:
        std.append( df_std[name + '_' + key][0])
    print(name)

    axs.bar(x + i * bwidth - bwidth * (n_cases - 1) / 2, mean, bwidth, yerr=std, label=labels[i], color=colors[i], alpha=0.7)

# Add labels, title, and custom ticks
axs.set_xlabel('Rate type')
axs.set_ylabel(r'Rate $s^{-1}$')
axs.set_title('TORC Plasmid')
axs.set_xticks(x)
axs.set_xticklabels(promoter_keys, rotation=45, ha='right')
axs.legend(loc='best')
axs.grid(True)

#plt.savefig(out_file+'png')
plt.show()
