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
dir_source = '../optimization/'
#v_code = 'test-dist_op_TORC_plasmid'
v_code = 'block-dist_op_TORC_plasmid_v2-01'

trials_file = dir_source + v_code + '-trials.pkl'
pkl_file  = dir_source + v_code+ '.pkl'
params_file = dir_source + v_code + '.csv'
loss_file = dir_source + v_code + '-values.csv'

out_file = 'loss_'+v_code

percentage_threshold = .10
#err_threshold = 0.5 # The minimum error we want


# Load
#-----------------------------------------------------------------------------------------------------------------------
with open(trials_file, 'rb') as file:
    trials_data = pickle.load(file)

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

# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we load
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True)


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

plt.show()
