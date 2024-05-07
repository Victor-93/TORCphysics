from TORCphysics import Circuit
import matplotlib.pyplot as plt
import seaborn as sns
# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# We want to run a single simulation of the Pleu500 with stochastic topoisomerase activities

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = None
enzymes_filename = None
environment_filename = 'environment_test.csv'
# environment_filename = 'environment.csv'
output_prefix = 'noRNAP'
frames = 1000#5000 #50000
series = True
continuation = False
dt = 2#0.5 #  0.25

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)
my_circuit.print_general_information()
#my_circuit.run()
enzymes_df, sites_df, environmental_df = my_circuit.run_return_dfs()

my_circuit.print_general_information()

# Plot stuff
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
nbins = [83,20]

fig, axs = plt.subplots(2, figsize=(width, 2*height), tight_layout=True)
for p, name in enumerate(names):

    ax = axs[p]

    print(name)

    # Load
    mask = enzymes_df['name'] == name
    x = enzymes_df[mask]['position'].to_numpy()

    # Plot
    # hist = sns.histplot(x, kde=True,  ax=ax, color=colors_dict[name], label=name)
    hist = sns.histplot(x, kde=True, bins=nbins[p], ax=ax, color=colors_dict[name], label=name)


    # Labels
    # ------------------------
    ax.set_ylabel('Density', fontsize=font_size)
    ax.set_xlabel(r'Position (bp)', fontsize=font_size)
    ax.set_title(name, fontsize=title_size)


plt.show()
x=2
