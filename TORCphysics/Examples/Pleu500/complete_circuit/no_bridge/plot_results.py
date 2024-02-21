import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from TORCphysics import Circuit
from TORCphysics import visualization as vs

# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by run_simulation.

# Inputs
# ---------------------------------------------------------
# Let's initialize circuit
circuit_filename = '../circuit.csv'
sites_filename = '../sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'nobridge'
frames = 50000
series = True
continuation = False
dt = 0.25
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

# Plotting inputs
# ---------------------------------------------------------
fa = int(frames / 2)
fb = frames
bins = 20
# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'Raspberry': 'red', 'lac1': 'green',
               'lac2': 'green'}
kwargs = {'linewidth': 2, 'ls': '-'}

# Let's plot
# ---------------------------------------------------------
csites_df = my_circuit.name + '_' + output_prefix + '_sites_df.csv'
cenzymes_df = my_circuit.name + '_' + output_prefix + '_enzymes_df.csv'

sites_df = pd.read_csv(csites_df, sep=',')
enzymes_df = pd.read_csv(cenzymes_df, sep=',')
enzymes_df = pd.read_csv(my_circuit.name + '_' + output_prefix + '_enzymes_df.csv')

# Plot responses
# -------------------------------------
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
fig.suptitle(my_circuit.name + ' rates')
vs.plot_site_response_curves(my_circuit, axs=axs, colors=colors_dict, site_type='gene')
plt.savefig('responses.png')


# Plot steady state curves
# -------------------------------------
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
fig.suptitle(my_circuit.name + ' rates')
vs.plot_steady_state_initiation_curve(my_circuit, sites_df, axs, fa=20000, fb=frames,
                                      colors=colors_dict, site_type='gene')  # fa=frames*.5, fb=frames,
plt.savefig('rates.png')

# Plot distribution of supercoiling
# -------------------------------------

# A bit of processing
mask = sites_df['type'] == 'gene'
names1 = sites_df[mask].drop_duplicates(subset='name')['name']
mask = sites_df['type'] == 'lacOP'
names2 = sites_df[mask].drop_duplicates(subset='name')['name']
names = pd.concat([names1, names2]).tolist()
time = np.arange(0, my_circuit.dt * (my_circuit.frames + 1), my_circuit.dt)

# And plot the superhelical density at sites
fig, axs = plt.subplots(len(names) + 1, figsize=(width, height * len(names) + 1), tight_layout=True)

# First for the global superhelical density
ax = axs[0]
name = 'circuit'
mask = sites_df['type'] == name
superhelical = sites_df[mask]['superhelical'].to_numpy()
sns.histplot(superhelical, kde=True, bins=50, ax=ax, color='black', label=name)

ax.set_ylabel('Density', fontsize=15)
ax.set_xlabel(r'Supercoiling density $(\sigma)$', fontsize=15)
ax.set_title('Global supercoiling', fontsize=15)

# Then all the other sites
for i, name in enumerate(names):
    ax = axs[i + 1]
    mask = sites_df['name'] == name
    superhelical = sites_df[mask]['superhelical'].to_numpy()
    sns.histplot(superhelical, kde=True, bins=30, ax=ax, color=colors_dict[name], label=name)

    ax.set_ylabel('Density', fontsize=15)
    ax.set_xlabel(r'Supercoiling density $(\sigma)$', fontsize=15)
    ax.set_title(name, fontsize=15)
#    axs.plot(time, superhelical, color=colors_dict[name], label=name, alpha=0.5, **kwargs)
plt.savefig('supercoiling_distribution.png')
