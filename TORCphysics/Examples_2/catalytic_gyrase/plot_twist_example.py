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
file_out = 'twist_example'
circuit_length = 2200

n_sim_shown = 4 # How many cases to plot

dt = 1.0 #0.25
initial_time = 0
final_time = 10800 // 10
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
    ax2 = ax.twinx()  # Different scale for RNAP
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
        # Linking dif
        ax = axs[n, 1]
        ax.plot(x, dLk[1:], lw=lw, color=colors[n_sim])
        ax2.plot(x, sigma[1:], lw=lw, color=colors[n_sim])
        # ax.grid(True)

    # Collect twists induced per binding event
    mtwists = []
    for n_sim in range(n_sims):

        # Linkding difference part and all the derivatives. These are topology variables
        sites_df = e_dict['sites_df'][n_sim]
        mask = sites_df['type'] == 'circuit'
        sigma = sites_df[mask]['superhelical'].to_numpy()
        dLk = sigma*Lk0 # This is the linking difference
        delta_dlK = np.diff(dLk[1:])/dt
        dif = np.diff(dLk[1:])
        dy = delta_dlK

        # And counting the number of enzymes - In case we need it
        sites_df = e_dict['sites_df'][n_sim]
        mask = sites_df['name'] == 'DNA_' + e_name
        enzyme = sites_df[mask].drop_duplicates('frame')
        nenzyme = enzyme['#enzymes'].to_numpy()#

        # Accumulate duration
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


    ax = axs[n, 0]
    ax.grid(True)
    #ax.set_ylabel(r'Density', fontsize=xlabel_size)
    ax.set_ylabel(r'Probability', fontsize=xlabel_size)
    ax.set_title(titles[n], fontsize=title_size)
    #if n ==0: # Topoisomerase
    #    ax.set_xlim(0,6)
    #else:
    #    ax.set_xlim(-8,0)

    # Histogram
    #ax.hist(mtwists, bins=nbins, density=True, color='blue', edgecolor='black')#, nbins=10)
    #sns.histplot(mtwists, kde=True, bins=nbins, ax=ax, color=hist_color, stat='probability')
    #sns.histplot(mtwists, kde=True, ax=ax, color=hist_color, stat='probability')
    sns.histplot(mtwists, kde=True, ax=ax, bins=10, color=hist_color, stat='probability')
axs[2,1].set_xlabel('Time (mins)', fontsize=xlabel_size)
axs[2,0].set_xlabel('Rotations induced per binding event', fontsize=xlabel_size)

#plt.savefig(file_out+'.png')
#plt.savefig(file_out+'.pdf')

plt.show()