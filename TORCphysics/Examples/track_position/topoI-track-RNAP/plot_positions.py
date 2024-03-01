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
out = 'TopoIRNAPTracking'
bins = 50
# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4

names = ['RNAP', 'topoI', 'gyrase']
colors_dict = {'RNAP': 'black', 'topoI': 'red', 'gyrase': 'cyan'}
kwargs = {'linewidth': 2, 'ls': '-'}
nbins = [40,80,80]

# Let's plot
# ---------------------------------------------------------
for p, name in enumerate(names):

    # Load
    x = np.loadtxt('position_'+name+'.txt')

    # Plot
    fig, ax = plt.subplots(1, figsize=(width, height), tight_layout=True)
    sns.histplot(x, kde=True, bins=nbins[p], ax=ax, color=colors_dict[name], label=name)
    ax.set_ylabel('Density', fontsize=15)
    ax.set_xlabel(r'Position (bp)', fontsize=15)
    ax.set_title(name, fontsize=15)

    plt.savefig('z_'+name+'.png')


