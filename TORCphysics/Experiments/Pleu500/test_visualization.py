from TORCphysics import Circuit
# from TORCphysics.src import analysis as an
import pandas as pd
# import numpy as np
import matplotlib.pyplot as plt
from TORCphysics import visualization as vs

# Description
# ---------------------------------------------------------
# This script is produced with the intention to test the visualization module.
# Even though this is not a formal test, we use it to check that we are producing the desired visualizations.

# Inputs
# ---------------------------------------------------------
csites_df = 'Pleu500_0_sites_df.csv'
cenzymes_df = 'Pleu500_0_enzymes_df.csv'
cenvironment_df = 'Pleu500_0_environment_df.csv'

log_file = 'Pleu500_0.log'

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

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'Raspberry': 'red'}
kwargs = {'linewidth': 2, 'ls': '--'}

# Load inputs
# ---------------------------------------------------------
sites_df = pd.read_csv(csites_df, sep=',')

# Create Figure
# ---------------------------------------------------------
fig, axs = plt.subplots(6, figsize=(width, 6 * height), tight_layout=True)

# Plot site response curves
# ---------------------------------------------------------
ax = axs[0]
vs.plot_site_response_curves(my_circuit, ax, site_type='gene', colors=colors_dict, **kwargs)

# Topoisomerase activity curves
# ---------------------------------------------------------
ax = axs[1]
vs.plot_topoisomerase_activity_curves_continuum(my_circuit, ax)

# Plot signal profiles
# ---------------------------------------------------------
ax = axs[2]
vs.plot_signal_profiles(my_circuit, sites_df, ax, site_type='gene', colors=colors_dict)

# Plot superhelical profiles
# ---------------------------------------------------------
ax = axs[3]
vs.plot_supercoiling_profiles(my_circuit, sites_df, ax, colors=colors_dict, site_type='gene')

# Plot cross-correlations
# ---------------------------------------------------------
ax = axs[4]
vs.plot_cross_correlation_with_site(my_circuit, sites_df, 'tetA', ax, fa=6000, fb=-1,
                                    ignore=['CDS'], site_type='gene', colors=colors_dict)

# Plot steady state initiation curve
# ---------------------------------------------------------
ax = axs[5]
vs.plot_steady_state_initiation_curve(my_circuit, sites_df, ax, ignore='CDS',
                                      fa=6000, fb=-1, colors=colors_dict, site_type='gene')

plt.savefig('visualization_test.png')
