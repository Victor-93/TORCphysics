from TORCphysics import Circuit
import pandas as pd
import matplotlib.pyplot as plt
from TORCphysics import visualization as vs
import numpy as np

# Initialise figure and circuit initial conditions
# ----------------------------------------------------------------------------------------------------------------------

# Figure
width = 8
height = 3
figure_output = 'calibration.png'

fig, axs = plt.subplots(2, figsize=(width, 2.5 * height), tight_layout=True)

# Circuit conditions
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_continuum_filename = 'environment_continuum.csv'
environment_filename = 'environment.csv'
output_prefix = 'continuum'
frames = 1000
series = True
continuation = False
mm = 'uniform'
dt = 1.0

csites_continuum_df = 'topos_continuum_sites_df.csv'
csites_stochastic_df = 'topos_stochastic_sites_df.csv'
# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------
continuum_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_continuum_filename,
                            output_prefix, frames, series, continuation, dt, 'continuum', mm)
stochastic_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, 'stoschastic', mm)

sites_continuum_df = pd.read_csv(csites_continuum_df, sep=',')
sites_stochastic_df = pd.read_csv(csites_stochastic_df, sep=',')

# Plot global supercoiling responses
# ---------------------------------------------------------
time = np.arange(0, (frames + 1) * dt, dt)

ax = axs[0]
mask = sites_continuum_df['type'] == 'circuit'  # This one contains global superhelical density
sigma_continuum = sites_continuum_df[mask]['superhelical'].to_numpy()
mask = sites_stochastic_df['type'] == 'circuit'  # This one contains global superhelical density
sigma_stochastic = sites_stochastic_df[mask]['superhelical'].to_numpy()
ax.plot(time, sigma_continuum, 'black')
ax.plot(time, sigma_stochastic, 'red')

#vs.plot_supercoiling_profiles(continuum_circuit, sites_continuum_df, ax, only_global=True, {'color':'black'})
#vs.plot_supercoiling_profiles(stochastic_circuit, sites_stochastic_df, ax, only_global=True, k{'color':'red'})
# Plot site response curves
# ---------------------------------------------------------
ax = axs[1]
vs.plot_site_response_curves(stochastic_circuit, ax)

plt.savefig('calibration.png')
