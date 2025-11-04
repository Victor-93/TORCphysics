import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import pickle
from TORCphysics import visualization as vs
from TORCphysics import analysis as an
# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# Let's plot the data from the example by loading the dictionaries with pickle

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
input_data = 'genearch_example.pkl'
file_out = 'genearch_example'
promoter_case = 'weak'

n_sim_shown = 4 # How many cases to plot

dt=1.0
initial_time = 0
final_time = 5400 #~1.5hrs
time = np.arange(initial_time, final_time + 2*dt, dt)
# -----------------------------------
# FIGURE Params
# -----------------------------------
width = 7
height = 3.5
lw = 3
font_size = 12
xlabel_size = 16
title_size = 18
nbins=30

colors = ['black', 'red', 'blue', 'green', 'yellow', 'purple', 'gray']

# ----------------------------------------------------------------------------------------------------------------------
# Load simulations
# ----------------------------------------------------------------------------------------------------------------------
with open(input_data, 'rb') as file:
    # Load the pickled data
    data = pickle.load(file)

# Unpack dictionary
sites_df_list = data['sites_df']
enzymes_df_list = data['enzymes_df']
environment_df_list = data['environment_df']

n_sims = len(sites_df_list) # How many simulations we launched
# ----------------------------------------------------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(2, figsize=(width, 2*height), tight_layout=True)#, sharex=True)

x = time / 60.0
for n in range(n_sim_shown): # Go through simulations

    # Select dataframes
    sites_df = sites_df_list[n]
    enzymes_df = enzymes_df_list[n]
    environment_df = environment_df_list[n]

    # Plot global supercoiling
    # ---------------------------------
    ax = axs[0]

    # Get global supercoiling
    mask = sites_df['type'] == 'circuit' # Global supercoiling can be obtained by type
    global_supercoiling = sites_df[mask]['superhelical']

    # Get local supercoiling
    mask = sites_df['name'] == promoter_case  # Local supercoiling can be obtained by site's name
    local_supercoiling = sites_df[mask]['superhelical']

    # And plot
    ax.plot(x, global_supercoiling, '-', lw=1, color=colors[n])
    ax.plot(x, local_supercoiling, '--', lw=1, color=colors[n])

    # Plot mRNA
    # ---------------------------------

    # Let's plot the transcription rate. We can obtain it by dividing the number of mRNA produced by the time
    gene_mask = environment_df['name'] == promoter_case # This gives the positions for which the condition is met
    mRNA = environment_df.loc[gene_mask] # This is the filtered dataframe

    axs[1].plot(mRNA['time']/60., mRNA['concentration']/mRNA['time']) # Plot the number of accumulated mRNAs/time (mRNA should degrade in real life, but we don't include it in the model)

    average_rate = np.mean(mRNA['concentration']/mRNA['time']) # Very vaguely, we can calculate the transcription rate as an average. Ideally we should avoid the first spikes....
    average_rate_s = "{:.2e}".format(average_rate)  # Let's save it as a string

    print('simulation ', n, 'transcription rate ', average_rate_s)

# Sort labels
axs[0].set_ylabel('Superhelical density')
axs[0].grid(True)

axs[1].set_ylabel(r'mRNA/time (second$^{-1}$)')
axs[1].grid(True)
axs[1].set_xlabel('Time (mins)', fontsize=xlabel_size)

#plt.savefig(file_out+'.png')
#plt.savefig(file_out+'.pdf')

plt.show()