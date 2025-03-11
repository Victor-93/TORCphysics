import numpy as np
import matplotlib.pyplot as plt


# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by the parallelization.

# Inputs
# ---------------------------------------------------------
out = 'TopoIRNAPTracking'

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

model_names = ['track-uniform', 'notrack-uniform']
colors_dict = {'track-uniform': 'black', 'notrack-uniform': 'red'}
kwargs = {'linewidth': 2, 'ls': '-'}
enzyme_names = ['RNAP', 'topoI', 'gyrase']
outside_label = ['a)', 'b)', 'c)']

# Let's plot
# ---------------------------------------------------------
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True)

# TODO: For the moment, for RNAP let's just plot the density, not the enrichment
for i, ename in enumerate(enzyme_names):
    ax = axs[i]

    # Get reference (topo activity when gene is turned off)
    if ename != 'RNAP':
        cdata = 'noRNAP/histogram_' + ename + '.txt'
        data = np.loadtxt(cdata)
        s0 = 0
        ref_y = data[s0:, 1]   # This is the reference data
    else:
        s0 = 130


    for p, name in enumerate(model_names):

        # Load data
        cdata = name + '/histogram_' + ename + '.txt'
        data = np.loadtxt(cdata)
        x = data[s0:, 0]
        dat_y = data[s0:, 1]   # This is our data

        # Transform data: This says the enrichment, hopefully?
        if ename != 'RNAP':
            y = dat_y / ref_y
        else:
            y = dat_y
        # print(ename, len(y))

        if ename == 'RNAP':
            RNAP_x = x
            RNAP_y = y
        if ename == 'topoI':
            topo_x = x
            topo_y = y

        # Plot results
        ax.plot(x, y, lw=lw, color=colors_dict[name], label=name)

    # Labels
    # ------------------------
    ax.set_ylabel('Enrichment', fontsize=font_size)
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    ax.set_title(ename, fontsize=title_size)
    ax.grid(True)
    ax.set_xlim(0, 10000)
    if ename != 'RNAP':
        ax.set_ylim(0, 2)

    # Add label outside the plot
    ax.text(-0.12, 0.99, outside_label[i], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')


axs[0].legend(loc='best', fontsize=font_size)

#plt.savefig(out+'.png')
#plt.savefig(out+'.pdf')
plt.show()

