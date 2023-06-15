from TORCphysics import Circuit
from TORCphysics.src import analysis as an
from TORCphysics.src import visualization as vs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Inputs
# ---------------------------------------------------------
csites_df = 'Pleu500_out_sites_df.csv'
cenzymes_df = 'Pleu500_out_enzymes_df.csv'
cenvironment_df = 'Pleu500_out_environment_df.csv'
sites_df = pd.read_csv(csites_df, sep=',')
enzymes_df = pd.read_csv(cenzymes_df, sep=',')

log_file = 'Pleu500_out.log'

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'Raspberry': 'red'}

circuit_filename = '../../../circuit.csv'
sites_filename = 'sites_maxmin.csv'
#sites_filename = 'sites_sam.csv'
enzymes_filename = '../../../enzymes.csv'
environment_filename = 'environment_stochastic.csv'
output_prefix = 'out'
frames = 500
series = True
continuation = False
tm = 'stochastic'
mm = 'uniform'
dt = .5

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt, tm, mm)
# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4

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
fig, axs = plt.subplots(1, figsize=(width,  height), tight_layout=True)

# Sites rate curves - Let's plot the rates modulated by supercoiling
# ---------------------------------------------------------
ax = axs
vs.plot_site_response_curves(my_circuit, ax)
ax.set_title('Response curves')
plt.savefig('responses.png')

# Plot cross-correlations
# ---------------------------------------------------------
fig, axs = plt.subplots(1, figsize=(width,  height), tight_layout=True)
ax = axs
vs.plot_signal_profiles(my_circuit, sites_df, ax, site_type='gene', colors=colors_dict)
ax.set_title('Signal profiles - OFF')
plt.savefig('signals.png')
