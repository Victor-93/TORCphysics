import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by the parallelization.

# Inputs
# ---------------------------------------------------------
circuit_name = 'lineartest'
n_sims = 10
out = 'noTracking'
bins = 50
# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

names = ['RNAP', 'topoI', 'gyrase']
colors_dict = {'RNAP': 'black', 'topoI': 'red', 'gyrase': 'cyan'}
kwargs = {'linewidth': 2, 'ls': '-'}
#nbins = [40,82,82]
nbins = [40,166,62]

# Let's plot
# ---------------------------------------------------------
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True)
for p, name in enumerate(names):

    ax = axs[p]

    # Load
    x = np.loadtxt('position_'+name+'.txt')

    # Plot
    sns.histplot(x, kde=True, bins=nbins[p], ax=ax, color=colors_dict[name], label=name)
    ax.set_ylabel('Density', fontsize=font_size)
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    ax.set_title(name, fontsize=title_size)
    ax.set_xlim(0, 10000)

plt.savefig(out+'.png')
plt.savefig(out+'.pdf')


