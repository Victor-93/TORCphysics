import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import pickle
import objective_single_molecule_experiment as osme

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This script is intended to analyse the results from the calibration process

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
#input_data = 'gyrase_calibration-trials.pkl'
input_data = 'gyrase_re-run.pkl'
file_out = 'twist_example'

# Forces tested during the calibration
Forces = [0.35, 0.6, 0.8, 1.0, 1.3] # pN

# -----------------------------------
# FIGURE Params
# -----------------------------------
width = 7
height = 4# 3.5
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
    best_case = pickle.load(file) # This file contains data already processed
# ----------------------------------------------------------------------------------------------------------------------
# Plot: Let's plot the distribution of induced rotations and the expected number of cycles.
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(2*width, height), tight_layout=True)#, sharex=True)

# Plot distribution of induced rotations per binding event for the best case
# -------------------------------------------------
# We need to flatten the list  (one entry per force, but we want it mixed)
rotations_dist = np.array([x for xs in best_case['rotations_distribution'] for x in xs])

# Play with the filter and ranges
rotations_dist = rotations_dist[(rotations_dist >= -20) & (rotations_dist <= 0)]
#sns.histplot(rotations_dist, kde=True, ax=axs[0], bins=21, color=hist_color, stat='probability')
sns.histplot(rotations_dist, kde=True, ax=axs[0], bins=20, binrange=(-20, 0), color=hist_color, stat='probability')

# Plot the expected number of cycles as a function of force
# -------------------------------------------------
x_force = np.arange(0.2,1.3,0.01) # continues force to calculate the theoretical curve
p_cycles, n_cycles = osme.mechano_gyrase(x_force) # This are the theoretical curves from bustamante's experiment
axs[1].plot(x_force,n_cycles)

# Calculated expected number of cycles as a function of force from simulation data
for i, rotations in enumerate(best_case['rotations_distribution']):
    n_rotations = len(rotations)
    n_cycles_sim = -np.mean(np.array(rotations) / 2.0)
    n_cycles_std = np.std(np.array(rotations) / 2.0) # standard deviation
    axs[1].errorbar(Forces[i], n_cycles_sim, yerr=n_cycles_std, fmt='-o', color='red')
    # axs[1].plot(Forces[i], n_cycles_sim, 'o', color='red')

axs[0].set_ylabel('Probability', fontsize=xlabel_size)
axs[0].set_xlabel('Induced rotations per binding event', fontsize=xlabel_size)

axs[1].set_xlabel('Force (pN)', fontsize=xlabel_size)
axs[1].set_ylabel('Expected number of cycles', fontsize=xlabel_size)


#plt.savefig(file_out+'.png')
#plt.savefig(file_out+'.pdf')

plt.show()