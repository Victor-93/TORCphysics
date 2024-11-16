import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by the parallelization.

# Inputs
# ---------------------------------------------------------
out = 'noTracking'
dt=1.0
# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

names = ['topoI', 'gyrase']
#names = ['RNAP', 'topoI', 'gyrase']
colors_dict = {'topoI': 'red', 'gyrase': 'cyan'}
#colors_dict = {'RNAP': 'black', 'topoI': 'red', 'gyrase': 'cyan'}
kwargs = {'linewidth': 2, 'ls': '-'}
#nbins = [166,62]
nbins = [301,201]

#tarray = [np.arange(0,6001,20), np.arange(0,6001,30)]
# Let's plot
# ---------------------------------------------------------
fig, axs = plt.subplots(2, figsize=(width, 2*height), tight_layout=True)
for p, name in enumerate(names):

    ax = axs[p]

    print(name)

    # Load
    #x = np.loadtxt('position_'+name+'.txt')
    x = np.loadtxt('position_'+name+'_dt'+str(dt)+'.txt')

    # Plot
    #hist = sns.histplot(x, ax=ax, color=colors_dict[name], label=name)

    hist = sns.histplot(x, kde=True, bins=nbins[p], ax=ax, color=colors_dict[name], label=name)

    # Extract info
    # ------------------------
    # Extract histogram data
    hist_data = hist.get_lines()[0].get_data()

    # Save histogram data
    hist_x, hist_y = hist_data
    #ax.plot(hist_x, hist_y, color='black', label='hist')

    # Save histogram data
    histogram_data = np.column_stack((hist_x, hist_y))
    np.savetxt('histogram_' + name + '.txt', histogram_data)

    # Plot KDE separately
    #kdeplot = sns.kdeplot(x)

    # Extract KDE data
    #kde_data = kdeplot.get_lines()[0].get_data()

    # Save KDE data
    #kde_x, kde_y = kde_data
    #ax.plot(kde_x, kde_y, color='black', label='dens')

    # Save KDE data
    #kde_curve_data = np.column_stack((kde_x, kde_y))
    #np.savetxt('kde_' + name +'.txt', kde_curve_data)

    # Labels
    # ------------------------
    ax.set_ylabel('Density', fontsize=font_size)
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    ax.set_title(name, fontsize=title_size)
    #ax.set_ylim(0,np.max(kde_y))
    #ax.set_xlim(0, 10000)

plt.show()
#plt.savefig(out+'.png')
#plt.savefig(out+'.pdf')


