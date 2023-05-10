from TORCphysics import Circuit
import matplotlib.pyplot as plt
from TORCphysics import visualization as vs
import sys

# TODO:
#  1.- Make the code work with stochastic topo binding. DONE
#  2.- Run experiment to test that it works. DONE
#  3.- Calibrate it so it works with Sam Meyers model, for now...
#  3.1.- I need to plot the global supercoiling response of the continuum case
#  3.2.- Plot the response of my initial case
#  3.3.- Run multiple experiments using the curve fit to find the optimal parameters of my stochastic topo binding

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# I will be plotting as the experiment advances.

# Initialise figure and circuit initial conditions
# ----------------------------------------------------------------------------------------------------------------------

# Figure
width = 8
height = 3
figure_output = 'calibration.png'

fig, axs = plt.subplots(3, figsize=(width, 3 * height), tight_layout=True)

# Circuit conditions
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment_continuum.csv'
output_prefix = 'continuum'
frames = 1500
series = True
continuation = False
mm = 'uniform'
dt = 1.0

# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------
tm = 'continuum'
continuum_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                            output_prefix, frames, series, continuation, dt, tm, mm)
continuum_circuit.run()

# Plot global supercoiling responses
ax = axs[0]
vs.plot_supercoiling_profiles(continuum_circuit, continuum_circuit.sites_df, ax, only_global=True)

# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------
tm = 'stochastic'
output_prefix = 'stochastic'
environment_filename = 'environment.csv'
n_simulations = 1
for ns in range(n_simulations):

    # Load simulation
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt, tm, mm)

    my_circuit.name = my_circuit.name + '_' + str(ns)
    my_circuit.log.name = my_circuit.name
    #  my_circuit.print_general_information()
    my_circuit.run()


# Plot site response curves
# ---------------------------------------------------------
ax = axs[1]
vs.plot_site_response_curves(my_circuit, ax)  #  , ignore=to_ignore)
# Let's ignore the

# Plot global supercoiling responses
# ---------------------------------------------------------
ax = axs[2]
vs.plot_supercoiling_profiles(my_circuit, my_circuit.sites_df, ax, only_global=True)

plt.savefig('calibration.png')
