import numpy as np
import matplotlib.pyplot as plt

# Description
# ---------------------------------------------------------
# I will process the densities...

# Inputs
# ---------------------------------------------------------
out = 'noTracking'

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

names = ['topoI', 'gyrase']
colors_dict = {'topoI': 'red', 'gyrase': 'cyan'}
kwargs = {'linewidth': 2, 'ls': '-'}

# Let's plot
# ---------------------------------------------------------
fig, axs = plt.subplots(2, figsize=(width, 2 * height), tight_layout=True)

for i, name in enumerate(names):
    ax = axs[i]

    # data =np.loadtxt('kde_'+name+'.txt')
    data = np.loadtxt('histogram_' + name + '.txt')

    x = data[:, 0]
    y = data[:, 1]
    ax.plot(x, y)

    # Labels
    # ------------------------
    ax.set_ylabel('Density', fontsize=font_size)
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    ax.set_title(name, fontsize=title_size)


plt.savefig('densities_count.png')