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
n_sims = 48
out = ['uniform', 'stall']
paths = ['RNAP_uniform/', 'RNAP_stall/'] #['no_bridge/multiple_runs/', 'bridge/multiple_runs/']
dt = 0.25
frames = 2000
fa = int(frames / 2)
fb = frames
bins =20
# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'Raspberry': 'red'}
kwargs = {'linewidth': 2, 'ls': '-'}

# Let's plot
# ---------------------------------------------------------
for p, path in enumerate(paths):

    # Let's filter positions.
    topoI =[]
    RNAP = []
    gyrase = []
    for n in range(n_sims):
        enzymes_df = pd.read_csv(path + circuit_name + '_' + str(n) + '__enzymes_df.csv')

        # RNAP
        mask = enzymes_df['name'] == 'RNAP'
        RNAP.append(enzymes_df[mask]['position'].values)

        # topoI
        mask = enzymes_df['name'] == 'topoI'
        topoI.append(enzymes_df[mask]['position'].values)

        # gyrase
        mask = enzymes_df['name'] == 'gyrase'
        gyrase.append(enzymes_df[mask]['position'].values)


    RNAP = np.concatenate(RNAP)
    topoI = np.concatenate(topoI)
    gyrase = np.concatenate(gyrase)
    print(len(RNAP), len(topoI), len(gyrase))


    #RNAP
    fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
    sns.histplot(RNAP, kde=True, bins=20, ax=axs, color='black', label='RNAP')
    axs.set_ylabel('Density', fontsize=15)
    axs.set_xlabel('Position', fontsize=15)
    axs.set_xlim(0, 5000)
    axs.set_title('RNAP')
    plt.savefig(out[p]+'_RNAP.png')


    #topoI
    fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
    sns.histplot(topoI, kde=True, bins=50, ax=axs, color='red', label='topoI')
    axs.set_ylabel('Density', fontsize=15)
    axs.set_xlabel('Position', fontsize=15)
    axs.set_xlim(0, 5000)
    axs.set_title('topoI')
    plt.savefig(out[p]+'_topoI.png')

    #gyrase
    fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
    sns.histplot(gyrase, kde=True, bins=50, ax=axs, color='cyan', label='gyrase')
    axs.set_ylabel('Density', fontsize=15)
    axs.set_xlabel('Position', fontsize=15)
    axs.set_xlim(0, 5000)
    axs.set_title('gyrase')
    plt.savefig(out[p]+'_gyrase.png')


