import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import pickle

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# Let's plot the data from the example by loading the dictionaries with pickle

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
input_data = 'example_data.pkl'
input_data = 'example_data_cycles.pkl'
file_out = 'twist_example'
circuit_length = 2200

bining_mode = 'active' # Bins the data according the time topoisomerases were active (twisting) on the DNA
bining_mode = 'bound'  # Bins the data according the periods of time topoisomerases remained bound to the DNA

n_sim_shown = 4 # How many cases to plot

dt = 1.0 #0.25
initial_time = 0
final_time = 10800
#final_time = 10800 // 2
time = np.arange(initial_time, final_time + dt, dt)
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
hist_color = 'blue'

colors = ['black', 'red', 'blue', 'green', 'yellow', 'purple', 'gray']

# ----------------------------------------------------------------------------------------------------------------------
# Load simulations
# ----------------------------------------------------------------------------------------------------------------------
with open(input_data, 'rb') as file:
    # Load the pickled data
    data = pickle.load(file)

# Unpack dictionary
gyrase_dict = data['gyrase']
gyrase_pos_dict = data['gyrase_positive']
topoI_dict = data['topoI']

n_sims = len(gyrase_dict['sites_df']) # How many simulations we launched
# ----------------------------------------------------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(3,2, figsize=(2*width, 3*height), tight_layout=True)#, sharex=True)

Lk0 = circuit_length/10.5

titles = ['Topoisomerase I Acting on -SC DNA', 'Gyrase Acting on RX DNA', 'Gyrase Acting on +SC DNA']
e_names = ['topoI', 'gyrase', 'gyrase']
window_size = 5 # Five seconds, which is the averaged time topoisomerases remain bound

#ylims = [ [-.12, 0.0], [-.12, .0], [-.12, .12] ]
ylims = [ [-24, 0], [-24, 0], [-24, 24] ]

for n, e_dict in enumerate([topoI_dict, gyrase_dict, gyrase_pos_dict]):

    e_name = e_names[n]

    # Let's first plot the linking number differences and superhelical density
    ax = axs[n, 1]
    ax.grid(True)
    ax2 = ax.twinx()  # Different scale for Superhelical density
    ax.set_ylabel(r'Linking Differnce', fontsize=xlabel_size)
    ax2.set_ylabel(r'Superhelical Density', fontsize=xlabel_size)
    ax.set_title(titles[n], fontsize=title_size)
    #ax.set_ylim(ylims[n])

    x = time/60.0
    for n_sim in range(n_sim_shown):
        # Linkding difference part and all the derivatives
        sites_df = e_dict['sites_df'][n_sim]
        mask = sites_df['type'] == 'circuit'
        sigma = sites_df[mask]['superhelical'].to_numpy()
        dLk = sigma*Lk0 # This is the linking difference
        ax = axs[n, 1]
        ax.plot(x, dLk[1:], lw=lw, color=colors[n_sim])  # Plots linking difference
        ax2.plot(x, sigma[1:], lw=lw, color=colors[n_sim])  # Plots superhelical density

    # Collect twists induced per binding event
    mtwists = []
    for n_sim in range(n_sims):

        # Linkding difference part and all the derivatives. These are topology variables
        sites_df = e_dict['sites_df'][n_sim]
        mask = sites_df['type'] == 'circuit'
        nenzymes = sites_df[mask]['#enzymes'].to_numpy()  #
        sigma = sites_df[mask]['superhelical'].to_numpy()
        dLk = sigma*Lk0 # This is the linking difference
        dif = np.diff(dLk[1:])  # Differences

        # Now much twist was introduced per cycle
        if bining_mode == 'active':
            xlabel = 'Rotations induced per active event'
            nonzero = np.abs(dif) > 1e-12
            diff = np.diff(nonzero.astype(int))
            p_start = np.where(diff==1)[0]+1
            p_end = np.where(diff==-1)[0]+1
            if nonzero[0]:
                p_starts = np.r_[0, p_start]
            if nonzero[-1]:
                p_end = np.r_[p_end, len(dif)]
            for start, end in zip(p_start, p_end):
                twist = dif[start:end].sum() # Summing the duration
                mtwists.append(twist)

        if bining_mode == 'bound':
            xlabel = 'Rotations induced per binding event'
            bound_df = sites_df[mask].drop_duplicates('frame')

            diff = np.diff(nenzymes[1:])
            p_start = np.where(diff == 1)[0] + 1
            p_end = np.where(diff == -1)[0] + 1

            for start, end in zip(p_start, p_end):
                twist = dif[start:end].sum() # Summing the duration
                mtwists.append(twist)

            # Boolean: enzyme active?
            #bound = bound_df['#enzymes'] == 1

            # Identify where activity changes
            # e.g. [0,0,1,1,1,0,0] → [0,0,1,1,1,2,2]
            #run_id = (bound != bound.shift()).cumsum()

            # Group by run_id, but keep only the active ones
            #groups = sites_df[mask].groupby(run_id)

            #mtwists = []
            #for rid, group in groups:
            #    if group['#enzymes'].iloc[0] == 1:
            #        dif = group['superhelical'].diff().fillna(0)
            #        mtwists.append(dif.sum() * Lk0)

    ax = axs[n, 0]
    ax.grid(True)
    ax.set_ylabel(r'Probability', fontsize=xlabel_size)
    ax.set_title(titles[n], fontsize=title_size)

    # Histogram
    if e_name == 'topoI':
        a = 0
        b = 5
    else:
        a = -5
        b = 0
    sns.histplot(mtwists, kde=True, ax=ax, bins=10, binrange=(a,b), color=hist_color, stat='probability')
    ax.set_xlim(a,b)

axs[2,1].set_xlabel('Time (mins)', fontsize=xlabel_size)
axs[2,0].set_xlabel(xlabel, fontsize=xlabel_size)

#plt.savefig(file_out+'.png')
#plt.savefig(file_out+'.pdf')

plt.show()