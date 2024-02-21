import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by comp_feed.py.
names = ['circuit', 'tetA', 'Raspberry', "mKalama1", 'lac1', 'lac2']
n_simulations = 96

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'Raspberry': 'red', 'lac1': 'green',
               'lac2': 'green', 'circuit': 'black'}
kwargs = {'linewidth': 2, 'ls': '-'}

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
    sns.histplot(superhelical, kde=True, bins=50, ax=ax, color=colors_dict[name], label=name)
    ax.set_ylabel('Density', fontsize=15)
    ax.set_xlabel(r'Supercoiling density $(\sigma)$', fontsize=15)
    if name == 'circuit':
        ctitle = 'global'
    else:
        ctitle = name
    ax.set_title(ctitle, fontsize=15)
plt.savefig('comp_feed-supercoiling_distribution.png')
plt.savefig('comp_feed-supercoiling_distribution.pdf')
