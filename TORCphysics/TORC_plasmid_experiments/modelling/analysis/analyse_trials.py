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
v_code = 'test-dist_op_TORC_plasmid'

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


# The rest is plotting the general loss... Maybe you can divide it per system...

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
df = df.sort_values(by='loss', ascending=False)#, inplace=True)
n = len(df['loss'])
nconsidered = int(n*percentage_threshold)
err_threshold = df['loss'].iloc[-nconsidered]
print('Number of tests', n)
print('Considered', nconsidered)
print('For ', percentage_threshold*100, '%')
# Filter according error
filtered_df = df[df['loss'] <= err_threshold]

loss = df['loss'].to_numpy()
floss = filtered_df['loss'].to_numpy()

# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we load
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True, sharex=True)

ms=6
ax = axs
ax.set_title('loss for '+ v_code)
#ax.plot(df['test'], df['loss'], 'o', ms=ms, color='blue', label='all')
#ax.plot(filtered_df['test'], filtered_df['loss'], 'o', ms=ms, color='red', label='best')

x = np.arange(1, n+1, 1)
ax.plot(loss, 'o', ms=ms, color='blue')
ax.plot(x[n-nconsidered:], loss[n-nconsidered:], 'o', ms=ms, color='red')

ax.grid(True)
ax.set_xlabel('test')
ax.set_ylabel('loss')
#ax.set_yscale('log')

#plt.savefig(out_file+'.png')
plt.show()
