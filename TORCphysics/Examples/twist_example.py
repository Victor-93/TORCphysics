import numpy as np
from TORCphysics import Circuit
import matplotlib.pyplot as plt
import random
import sys
import seaborn as sns
# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# Let's try to test a plasmid of the same size as the one used in the study, with one binding site, and
#
# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0 #0.25
initial_time = 0
final_time = 10800
#final_time = 10800 // 10
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'ref1_c05_twist-contribution'
n_sims = 100  # This is the total
n_sim_shown = 4

# For the simulation
circuit_filename = '../../circuit.csv'
circuit_filename = 'circuit_linear.csv' # This is the one they used in the experiment
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
environment_topoI = 'environment_topoI.csv'
environment_gyrase = 'environment_gyrase.csv'


# Superhelical values (sigma) for each case
sigma_0_topo = -0.11#0.075  # Approximately -20 supercoils according the paper
sigma_0_gyrase_relaxed = 0.0  # We suppose this one.
sigma_0_gyrase_positive = 0.11  # We suppose this one.

output_prefix = 'test0'
series = True
continuation = False

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
# Simulations
# ----------------------------------------------------------------------------------------------------------------------

# For gyrase Relaxed ------------------------------------
environment_filename = environment_gyrase
sigma0 = sigma_0_gyrase_relaxed
enzymes_df_list = []
sites_df_list = []

#n_site = 50 # The site with k_on different of 0
my_type = 'DNA_gyrase'
for n in range(n_sims):
    # Initialize circuit with the initial conditions
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt)
    my_circuit.name = my_circuit.name + '_' + str(n) # We can change the name of the circuit
    my_circuit.seed = my_circuit.seed + n + random.randrange(sys.maxsize)
    my_circuit.rng = np.random.default_rng(my_circuit.seed)
    for enzyme in my_circuit.enzyme_list:
        enzyme.superhelical = sigma0
    my_circuit.update_twist()
    my_circuit.update_supercoiling()
    my_circuit.update_global_twist()
    my_circuit.update_global_superhelical()

    n_sites = len(my_circuit.site_list)
    n_site = n_sites // 2 # Placed at the middle
    # Make binding sites with k_on == 0 except for n_site
    k_on = my_circuit.site_list[1].k_on
    binding_model = my_circuit.site_list[1].binding_model
    binding_oparams = my_circuit.site_list[1].binding_model_oparams
    for s, site in enumerate(my_circuit.site_list):
        site.k_on = 0.0
        if site.site_type == my_type and s >0:
#        if 'DNA_' not in site.name:
            site.binding_model_oparams['k_on'] = 0.0
            site.binding_model = None
            site.binding_model_oparams = None
    my_circuit.site_list[n_site].k_on = k_on
    my_circuit.site_list[n_site].binding_model = binding_model
    my_circuit.site_list[n_site].binding_model_oparams = binding_oparams

    # Run simulations but storing dataframes on memory then adding them to the lists.
    enzymes_df, sites_df, environment_df = my_circuit.run_return_dfs() # Function run_return_dfs() returns the dataframes with the results of the simulation (it does not write CSV files).

    # Append dataframes to the lists.
    enzymes_df_list.append(enzymes_df)
    sites_df_list.append(sites_df)

gyrase_dict = {'sites_df': sites_df_list}

# For gyrase Positive ------------------------------------
environment_filename = environment_gyrase
sigma0 = sigma_0_gyrase_positive
enzymes_df_list = []
sites_df_list = []

#n_site = 50 # The site with k_on different of 0
my_type = 'DNA_gyrase'
for n in range(n_sims):
    # Initialize circuit with the initial conditions
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt)
    my_circuit.name = my_circuit.name + '_' + str(n) # We can change the name of the circuit
    my_circuit.seed = my_circuit.seed + n + random.randrange(sys.maxsize)
    my_circuit.rng = np.random.default_rng(my_circuit.seed)
    for enzyme in my_circuit.enzyme_list:
        enzyme.superhelical = sigma0
    my_circuit.update_twist()
    my_circuit.update_supercoiling()
    my_circuit.update_global_twist()
    my_circuit.update_global_superhelical()

    n_sites = len(my_circuit.site_list)
    n_site = n_sites // 2 # Placed at the middle
    # Make binding sites with k_on == 0 except for n_site
    k_on = my_circuit.site_list[1].k_on
    binding_model = my_circuit.site_list[1].binding_model
    binding_oparams = my_circuit.site_list[1].binding_model_oparams
    for s, site in enumerate(my_circuit.site_list):
        site.k_on = 0.0
        if site.site_type == my_type and s >0:
#        if 'DNA_' not in site.name:
            site.binding_model_oparams['k_on'] = 0.0
            site.binding_model = None
            site.binding_model_oparams = None
    my_circuit.site_list[n_site].k_on = k_on
    my_circuit.site_list[n_site].binding_model = binding_model
    my_circuit.site_list[n_site].binding_model_oparams = binding_oparams

    # Run simulations but storing dataframes on memory then adding them to the lists.
    enzymes_df, sites_df, environment_df = my_circuit.run_return_dfs() # Function run_return_dfs() returns the dataframes with the results of the simulation (it does not write CSV files).

    # Append dataframes to the lists.
    enzymes_df_list.append(enzymes_df)
    sites_df_list.append(sites_df)

gyrase_pos_dict = {'sites_df': sites_df_list}

# For topoI ----------------------------------
environment_filename = environment_topoI
sigma0 = sigma_0_topo
enzymes_df_list = []
sites_df_list = []
my_type = 'DNA_topoI'

for n in range(n_sims):
    # Initialize circuit with the initial conditions
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt)
    my_circuit.name = my_circuit.name + '_' + str(n) # We can change the name of the circuit
    my_circuit.seed = my_circuit.seed + n + random.randrange(sys.maxsize)
    my_circuit.rng = np.random.default_rng(my_circuit.seed)

    n_sites = len(my_circuit.site_list)
    n_site = n_sites // 2 # Placed at the middle

    for enzyme in my_circuit.enzyme_list:
        enzyme.superhelical = sigma0
    my_circuit.update_twist()
    my_circuit.update_supercoiling()
    my_circuit.update_global_twist()
    my_circuit.update_global_superhelical()

    # Make binding sites with k_on == 0 except for n_site
    k_on = my_circuit.site_list[1].k_on
    binding_model = my_circuit.site_list[1].binding_model
    binding_oparams = my_circuit.site_list[1].binding_model_oparams
    for s, site in enumerate(my_circuit.site_list):
        site.k_on = 0.0
        if site.site_type == my_type and s > 0:
            #        if 'DNA_' not in site.name:
            site.binding_model_oparams['k_on'] = 0.0
            site.binding_model = None
            site.binding_model_oparams = None
    my_circuit.site_list[n_site].k_on = k_on
    my_circuit.site_list[n_site].binding_model = binding_model
    my_circuit.site_list[n_site].binding_model_oparams = binding_oparams

    # Run simulations but storing dataframes on memory then adding them to the lists.
    enzymes_df, sites_df, environment_df = my_circuit.run_return_dfs() # Function run_return_dfs() returns the dataframes with the results of the simulation (it does not write CSV files).

    # Append dataframes to the lists.
    enzymes_df_list.append(enzymes_df)
    sites_df_list.append(sites_df)

topoI_dict = {'sites_df': sites_df_list}

# ----------------------------------------------------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(3,2, figsize=(2*width, 3*height), tight_layout=True)#, sharex=True)

Lk0 = my_circuit.size/10.5

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

        # Linkding difference part and all the derivatives
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
    sns.histplot(mtwists, kde=True, ax=ax, color=hist_color, stat='probability')
axs[2,1].set_xlabel('Time (mins)', fontsize=xlabel_size)
axs[2,0].set_xlabel('Rotations induced per binding event', fontsize=xlabel_size)

plt.savefig(file_out+'.png')
plt.savefig(file_out+'.pdf')

plt.show()