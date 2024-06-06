import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from TORCphysics import Circuit
from TORCphysics import visualization as vs


# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by comp_feed.py.
names = ['circuit', 'tetA', "mKalama1", 'mRaspberry', 'antitet', 'bla']
n_simulations = 12
with_topoI = False

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = '../enzymes.csv'
if with_topoI:
    environment_filename = '../environment.csv'
    output_prefix = 'WT_system'
    enzyme_names = ['RNAP', 'topoI', 'gyrase']
else:
    environment_filename = '../environment_notopoI.csv'
    output_prefix = 'notopoI_system'
    enzyme_names = ['RNAP', 'gyrase']
frames = 5000 #50000
series = True
continuation = False
dt = 0.25

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'mRaspberry': 'red', 'lac1': 'green',
               'lac2': 'green', 'antitet':'green', 'circuit':'black', 'RNAP': 'black', 'topoI': 'red', 'gyrase':'cyan',
               'bla': 'purple'}
kwargs = {'linewidth': 2, 'ls': '-'}


names_genes = ['tetA', 'mKalama1', 'antitet', 'mRaspberry', 'bla']

# Binding
# -------------------------------------
fig, axs = plt.subplots(len(names_genes), figsize=(width, height * len(names_genes)), tight_layout=True)

out_file = output_prefix + '-bind'
for i, name in enumerate(names_genes):

    # Extract data
    # ------------------------------------------------------------------------------------
    print(name)
    input_df = output_prefix + '-binding_' + name + '_df.csv'
    bind_df = pd.read_csv(input_df, sep=',')

    time = bind_df['time'].to_numpy().ravel()

    # TODO: Plotea el rate como le hacias en el otro, dividido por el tiempo.
    sum_simu_array = bind_df.filter(like='simu').sum(axis=1).tolist()#.to_numpy()

    sum_ini = np.zeros_like(sum_simu_array)
    for k in range(1, len(time)):
        sum_ini[k] = np.sum(sum_simu_array[0:k+1])


    a = sum_ini[-1]
    t2 = time[-1]
    print(a, t2, a/t2)
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

plt.savefig(out_file+'.png')


# Plot distribution of supercoiling
# -------------------------------------
fig, axs = plt.subplots(len(names), figsize=(width, height * len(names)), tight_layout=True)

out_file = output_prefix + '-supercoiling'
for i, name in enumerate(names):

    # Extract data
    # ------------------------------------------------------------------------------------
    print(name)
    input_df = output_prefix + '-superhelical_' + name + '_df.csv'
    superhelical_df = pd.read_csv(input_df, sep=',')

    superhelical = pd.concat([superhelical_df[col] for col in superhelical_df.columns if col.startswith('simu')], axis=1)
    superhelical = superhelical.to_numpy().ravel()

    # Plot
    # ------------------------------------------------------------------------------------
    ax = axs[i]
    sns.histplot(superhelical, kde=True, bins=400, ax=ax, color=colors_dict[name], label=name)
    ax.set_ylabel('Density', fontsize=15)
    ax.set_xlabel(r'Supercoiling density $(\sigma)$', fontsize=15)
    if name == 'circuit':
        ctitle = 'global'
    else:
        ctitle = name
    ax.set_title(ctitle, fontsize=15)
    ax.set_xlim(-0.15, 0.1)

plt.savefig(out_file+'.png')

# Plot distribution of enzymes
# -------------------------------------
fig, axs = plt.subplots(len(enzyme_names), figsize=(width, height * len(enzyme_names)), tight_layout=True)

out_file = output_prefix + '-enzymes'
for i, name in enumerate(enzyme_names):

    factor = .5
    # Get environmental object
    my_environmental = [environmental for environmental in my_circuit.environmental_list if environmental.name == name][
        0]
    nbins = int(factor * my_circuit.size / my_environmental.size)  # number of bins. - note that this only applies for topos

    # Extract data
    # ------------------------------------------------------------------------------------
    print(name)
    input_df = output_prefix +'-positions_' + name + '.txt'
    x = np.loadtxt(input_df)

    # Plot
    # ------------------------------------------------------------------------------------
    ax = axs[i]
    sns.histplot(x, kde=True, bins=nbins, ax=ax, color=colors_dict[name], label=name)
    ax.set_ylabel('Density', fontsize=15)
    ax.set_xlabel(r'Position (bp)', fontsize=15)
    if name == 'circuit':
        ctitle = 'global'
    else:
        ctitle = name
    ax.set_title(ctitle, fontsize=15)
    #ax.set_xlim(-0.15, 0.1)

plt.savefig(out_file+'.png')

# Plot responses
# -------------------------------------
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
fig.suptitle(my_circuit.name + ' rates')
vs.plot_site_response_curves(my_circuit, axs=axs, colors=colors_dict, site_type='gene')
plt.savefig('responses.png')