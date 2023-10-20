import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from TORCphysics import visualization as vs
from TORCphysics import analysis as an
import sys

# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by the parallelization.

# Inputs
# ---------------------------------------------------------
circuit_name = 'Pleu500'
n_sims = 96
title = ['ON', 'OFF']
paths = ['no_bridge/multiple_runs/', 'bridge/multiple_runs/']
dt = 0.5
frames = 8000
fa = int(frames / 2)
fb = frames
ref_name = 'tetA'

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'Raspberry': 'red'}
kwargs = {'linewidth': 2, 'ls': '-'}

# Let's plot
# ---------------------------------------------------------
for p, path in enumerate(paths):
    fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)

    # Let's compute the cross correlations and save them in a new array
    cross_list = []
    for n in range(1, n_sims + 1):
        sites_df = pd.read_csv(path + circuit_name + '_' + str(n) + '_out_sites_df.csv')
        signals, names = an.build_signal_by_type(sites_df, 'gene')
        signals_t0 = []
        for signal in signals:
            signals_t0.append(signal[fa:fb])
        cross, lag = an.cross_correlation_hmatrix(signals_t0, dt)
        cross_list.append(cross)

    # Let's calculate averages
    m = len(signals_t0)
    b = len(signal[fa:fb])
    avg_cross = np.zeros((m,m,b))
    std_cross = np.zeros((m,m,b))
    # AVG
    for cross in cross_list:
        for i in range(m):
            for j in range(m):
                avg_cross[i, j, :] = avg_cross[i, j, :] + cross[i, j, :]
    avg_cross = avg_cross / n_sims
    # STD
    for cross in cross_list:
        for i in range(m):
            for j in range(m):
                std_cross[i, j, :] = std_cross[i, j, :]+np.square(cross[i, j, :] - avg_cross[i, j, :])
    std_cross = np.sqrt(std_cross / n_sims)/np.sqrt(n_sims)

    for i, name in enumerate(names):
        if name == ref_name:
            index = i
        if name == 'mKalama1':
            mKalama_index = i
        if name == 'Raspberry':
            Raspberry_index = i

    # Let's plot the correlations
    # ---------------------------------------
    # We need to find the maximum correlation write it
    indices = [mKalama_index, Raspberry_index]

    for i in indices:
        maxlag = lag[np.argmax(avg_cross[index, i, :])]
        my_label = names[i] + f' lag={maxlag:.2f}s'
        axs.plot(lag, avg_cross[index, i, :], color=colors_dict[names[i]], label=my_label, **kwargs)
        axs.fill_between(lag, avg_cross[index, i, :]-std_cross[index, i, :],
                         avg_cross[index, i, :]+std_cross[index, i, :],
                         color=colors_dict[names[i]], alpha=0.25)
    axs.set_ylabel('Correlation', fontsize=15)
    axs.set_xlabel('Time lag (seconds)', fontsize=15)
    axs.legend(loc='best')
    axs.set_title('Cross-Correlation with tetA - ' + title[p], fontsize=20)
    axs.grid(True)
    axs.set_xlim(-150, 150)
    plt.savefig(title[p]+'_v2.png')



