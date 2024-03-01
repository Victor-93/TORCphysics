import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by comp_feed.py.
names = ['circuit', 'tetA', "mKalama1", 'lac1', 'lac2']
n_simulations = 96

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'Raspberry': 'red', 'lac1': 'green',
               'lac2': 'green', 'circuit': 'black'}
kwargs = {'linewidth': 2, 'ls': '-'}


names_genes = ['tetA', 'mKalama1']
# Binding
# -------------------------------------
fig, axs = plt.subplots(len(names_genes), figsize=(width, height * len(names_genes)), tight_layout=True)

for i, name in enumerate(names_genes):

    # Extract data
    # ------------------------------------------------------------------------------------
    print(name)
    input_df = 'binding_' + name + '_df.csv'
    bind_df = pd.read_csv(input_df, sep=',')

    time = bind_df['time'].to_numpy().ravel()

    # TODO: Plotea el rate como le hacias en el otro, dividido por el tiempo.
    sum_simu_array = bind_df.filter(like='simu').sum(axis=1).tolist()#.to_numpy()

    sum_ini = np.zeros_like(sum_simu_array)
    for k in range(1, len(time)):
        sum_ini[k] = np.sum(sum_simu_array[0:k])

    #curve = np.log(sum_ini / time)
    curve = sum_ini/(time*n_simulations)

    # Plot
    # ------------------------------------------------------------------------------------
    ax = axs[i]
    ax.plot(time, curve, color=colors_dict[name], **kwargs)
    ax.set_ylabel('rate', fontsize=15)
    ax.set_xlabel(r'time (s)', fontsize=15)
    ctitle = name
    ax.set_title(ctitle, fontsize=15)
    ax.grid(True)
    #ax.set_ylim(0,.012)

plt.savefig('comp_feed-bind.png')
plt.savefig('comp_feed-bind.pdf')


# Plot distribution of supercoiling
# -------------------------------------
fig, axs = plt.subplots(len(names), figsize=(width, height * len(names)), tight_layout=True)

for i, name in enumerate(names):

    # Extract data
    # ------------------------------------------------------------------------------------
    print(name)
    input_df = 'superhelical_' + name + '_df.csv'
    superhelical_df = pd.read_csv(input_df, sep=',')

    superhelical = pd.concat([superhelical_df[col] for col in superhelical_df.columns if col.startswith('simu')], axis=1)
    superhelical = superhelical.to_numpy().ravel()

    # Plot
    # ------------------------------------------------------------------------------------
    ax = axs[i]
    sns.histplot(superhelical, kde=True, bins=200, ax=ax, color=colors_dict[name], label=name)
    ax.set_ylabel('Density', fontsize=15)
    ax.set_xlabel(r'Supercoiling density $(\sigma)$', fontsize=15)
    if name == 'circuit':
        ctitle = 'global'
    else:
        ctitle = name
    ax.set_title(ctitle, fontsize=15)
    ax.set_xlim(-0.07, 0.01)


plt.savefig('comp_feed-supercoiling_distribution.png')
plt.savefig('comp_feed-supercoiling_distribution.pdf')

