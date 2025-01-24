import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns


# Description
#-----------------------------------------------------------------------------------------------------------------------
# We want to plot the expression rates from the TORC plasmid following the optmization process

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
dir_source = '../optimization/'
v_code = 'block-full-dist_op_TORC_plasmid_st3.2-01'
#v_code = 'block-full-trackingON-dist_op_TORC_plasmid_st3.2-01'
#v_code = 'block-op_TORC_plasmid'
#v_code = 'block-rep-TORC_plasmid'
#v_code = 'op_TORC_plasmid'
#v_code = 'block-dist_op_TORC_plasmid'
#v_code = 'batch-dist_op_TORC_plasmid-01'
#v_code = 'full-dist_op_TORC_plasmid_st3-02'
#v_code = 'full-trackingON-dist_op_TORC_plasmid_st3-02'


pkl_file  = dir_source + v_code+ '.pkl'

out_file = 'superhelical_'+v_code

n_simulations = 120#24 #180 - ask for 61 cores (Each would run 3 simulations per system approximately).
n_batches = 12 # The number of simulations n_simulations will be grouped into n_batches. For each batch, an average
               # will be calculated. In this way, we will be comparing averages with averages.

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 12#8
height = 5
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

colors=['yellow', 'blue', 'red', 'green', 'purple']
#gene_names = ['PleuWT', 'tetA', 'antitet', 'bla']#, 'lacI']

# Load
#-----------------------------------------------------------------------------------------------------------------------
with open(pkl_file, 'rb') as file:
    pkl_data = pickle.load(file)


# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we do the process
fig, axs = plt.subplots(1, figsize=(width, 1*height), tight_layout=True, sharex=True)

ax = axs
ax.set_title('Version '+ v_code)

my_dict = {}
for system in pkl_data:
    name = system['name']
#    g = []
#    frames = len(system['global_superhelical'][0,:]) # How ma

#    for n in range(n_simulations):
#        g.append(system['global_superhelical'][n,int(frames/2):]) # Only the last half
#    my_dict[name] = np.array(g).flatten()
    my_dict[name] = system['global_superhelical'].flatten()

global_superhelical_df = pd.DataFrame.from_dict(my_dict)

sns.violinplot(data=global_superhelical_df, ax=ax, inner="quart")#, cut=0, color=colors[i])
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_ylabel('Global superhelical density')
ax.grid(True)
ax.set_ylim(-.11, -.0)

plt.show()



