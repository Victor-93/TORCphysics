from TORCphysics import analysis as an
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Description
# ---------------------------------------------------------
# This is script is intended as an example to show how to process the outputs of TORCphysics.
# We will plot: 1.- the binding/transcription signals. 2.- Cross-correlations. 3.- Binding curves.
# 4.- Global supercoiling & supercoiling at site. 5.- topoisomerase activity curves (continuum). 6.- Rates.
# 7.- Supercoiling distributions at sites.

# Inputs
# ---------------------------------------------------------
csites_df = 'Pleu500_0_sites_df.csv'
cenzymes_df = 'Pleu500_0_enzymes_df.csv'
cenvironment_df = 'Pleu500_0_environment_df.csv'

sites_csv = 'Pleu500_0_sites_output.csv'
environment_csv = 'Pleu500_0_environment_output.csv'

log_file = 'Pleu500_0.log'

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 3

# Better use these for the colors of genes...
# Sort them according the input file...
colors = []
colors.append("yellow")
colors.append("green")
colors.append("blue")
colors.append("red")
colors.append("magenta")
colors.append("black")
colors.append("cyan")
colors.append("black")


# Functions that will be useful
# ---------------------------------------------------------
def ax_params(axis, xl, yl, grid, legend):
    axis.grid(grid)
    axis.set_ylabel(yl)
    axis.set_xlabel(xl)
    if legend:
        axis.legend(loc='best')


# Load inputs
# ---------------------------------------------------------
sites_df = pd.read_csv(csites_df, sep=',')
#enzymes_df = pd.read_csv(cenzymes_df, sep=',')
dt = 1.0  # This should be extracted from the log file
# Create Figure
# ---------------------------------------------------------
fig, axs = plt.subplots(3, figsize=(width, 3 * height), tight_layout=True)

# Signals
# ---------------------------------------------------------
ax = axs[0]

signals, names = an.build_signal_by_type(sites_df, 'gene')
time = np.arange(0, dt * len(signals[0]), dt)

for i, signal in enumerate(signals):
    ax.plot(time, signal, color=colors[i], label=names[i], alpha=0.5)
ax_params(axis=ax, xl='time (seconds)', yl='Transcription signal', grid=True, legend=False)

# Superhelical density
# ---------------------------------------------------------
ax = axs[1]

# Let's plot the global superhelical density
mask = sites_df['type'] == 'circuit'  # This one contains global superhelical density
global_sigma = sites_df[mask]['superhelical'].to_numpy()

# And plot the superhelical density at sites
mask = sites_df['type'] == 'gene'
genes_names = sites_df[mask].drop_duplicates(subset='name')['name']  # Let's filter the genes so we get the gene names
for i, name in enumerate(genes_names):
    mask = sites_df['name'] == name
    superhelical = sites_df[mask]['superhelical'].to_numpy()
    ax.plot(time, superhelical, color=colors[i], label=name)
ax.plot(time, global_sigma, color='black', label='global')
ax_params(axis=ax, xl='time (seconds)', yl='Supercoiling at site', grid=True, legend=True)

# Cross-correlations
# ---------------------------------------------------------
ax = axs[2]
t0 = 2000  # Time in which we assume the system has reached the steady state
# Signals from t0 to end
signals_t0 = []
for signal in signals:
    signals_t0.append(signal[t0:])
cross, lag = an.cross_correlation_hmatrix(signals_t0, dt)
for i, name in enumerate(genes_names):
    if name == 'tetA':
        continue
    ax.plot(lag, cross[0,i,:], color=colors[i], label=name)
ax_params(axis=ax, xl='time lag (seconds)', yl='Cross-correlation w tetA', grid=True, legend=True)

ax.set_xlim(-200, 200)

plt.savefig('signals.png')
