import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d


# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by the parallelization.

# Inputs
# ---------------------------------------------------------
out = 'track-uniform'
# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

names = ['RNAP', 'topoI', 'gyrase']
colors_dict = {'topoI': 'red', 'gyrase': 'cyan'}
colors_dict = {'RNAP': 'black', 'topoI': 'red', 'gyrase': 'cyan'}
kwargs = {'linewidth': 2, 'ls': '-'}
#nbins = [200, 166,62]
nbins = [40,166,62]


# Define ranges and x-axis for interpolation/plotting
gene_start = 4000
gene_end = 6000

nbp = 10000

x_spacing = 10

# Define x-axes
x_system = np.arange(1, nbp, x_spacing)
x_gene = np.arange(gene_start, gene_end, x_spacing)

# Let's plot
# ---------------------------------------------------------
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True)
for p, name in enumerate(names):

    ax = axs[p]

    # Load
    x = np.loadtxt('position_'+name+'.txt')

    x = x[~np.isnan(x)]  # Just in case to remove nans

    print(name, len(x))

    # Calculate the histogram without normalization
    counts, bin_edges = np.histogram(x, bins=nbins[p], density=False)

    # Calculate the bin width
    bin_width = bin_edges[1] - bin_edges[0]

    # calculate kde with scipy
    kde = gaussian_kde(x)
    kde_x = np.linspace(min(x), max(x), nbins[p])
    kde_y = kde(kde_x)

    # Scale KDE values
    kde_y_scaled = kde_y * len(x) * bin_width

    # Create interpolation function
    interp_fun = interp1d(kde_x, kde_y_scaled, kind='linear', fill_value='extrapolate')  # everything default

    if name == 'RNAP':
        x_interpolated = x_gene
    else:
        x_interpolated = x_system

    # Get interpolated y-values
    y_interpolated = interp_fun(x_interpolated)

    # Plot
    #hist = sns.histplot(x, ax=ax, color=colors_dict[name], label=name)

    hist = sns.histplot(x, kde=True, bins=nbins[p], ax=ax, color=colors_dict[name], label=name)
    ax.plot(kde_x, kde_y_scaled, color='green', linewidth=lw*1.5, linestyle='--') # Let's plot it to see if it changes
    ax.plot(x_interpolated, y_interpolated, color='purple', linewidth=lw*1.5, linestyle='--') # Let's plot it to see if it changes

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
    ax.set_xlim(0, 10000)

plt.savefig(out+'.png')
plt.savefig(out+'.pdf')


