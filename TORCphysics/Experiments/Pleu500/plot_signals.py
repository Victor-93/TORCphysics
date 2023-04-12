from TORCphysics import Circuit
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

log_file = 'Pleu500_0.log'

# circuit_csv = 'circuit.csv'
# sites_csv = 'Pleu500_0_sites_output.csv'
# enzymes_csv = 'Pleu500_0_enzymes_output.csv'
# environment_csv = 'Pleu500_0_environment_output.csv'
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'output'
frames = 10000
series = True
continuation = False
tm = 'continuum'
mm = 'uniform'
dt = 1.0
n_simulations = 1
bridge_time = 40000

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt, tm, mm)
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
# enzymes_df = pd.read_csv(cenzymes_df, sep=',')
dt = 1.0  # This should be extracted from the log file
# Create Figure
# ---------------------------------------------------------
fig, axs = plt.subplots(6, figsize=(width, 6 * height), tight_layout=True)

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
    ax.plot(lag, cross[0, i, :], color=colors[i], label=name)
ax_params(axis=ax, xl='time lag (seconds)', yl='Cross-correlation w tetA', grid=True, legend=True)

ax.set_xlim(-200, 200)

# Sites rate curves - Let's plot the rates modulated by supercoiling
# ---------------------------------------------------------
ax = axs[3]
i = -1
for site in my_circuit.site_list:
    if site.site_type == 'gene':
        i += 1
        rate, x = an.site_activity_curves(site)
        ax.plot(x, rate, color=colors[i], label=site.name)

ax_params(axis=ax, xl=r'\sigma', yl=r'Initiation rate ($s^{-1}$)', grid=True, legend=True)

# Topoisomerase activity curves - Let's plot the rates modulated by supercoiling
# ---------------------------------------------------------
ax = axs[4]
i = -1
for environmental in my_circuit.environmental_list:
    if environmental.enzyme_type == 'topo':
        i += 1

        topo_curve, x = an.topoisomerase_activity_curves_continuum(environmental, dt=dt)
        ax.plot(x, topo_curve, color=colors[i], label=environmental.name)
        if i == 0:
            topo_sum = np.zeros_like(topo_curve)
        topo_sum += topo_curve
ax.plot(x, topo_sum, color='black', label='sum')
ax_params(axis=ax, xl=r'\sigma', yl=r'$\sigma$ removed per timestep', grid=True, legend=True)

# Initiation rate
# ---------------------------------------------------------
ax = axs[5]
curves, rates, labels = an.initiation_rate_by_type(sites_df, 'gene', time, ta=5000, tb=-1)

for i, curve in enumerate(curves):
    ax.plot(time, curve, color=colors[i], label=names[i], alpha=0.5)
    print(rates[i])
ax_params(axis=ax, xl='time (seconds)', yl='curve', grid=True, legend=True)

plt.savefig('signals.png')
