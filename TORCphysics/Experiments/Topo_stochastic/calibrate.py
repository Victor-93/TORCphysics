from TORCphysics import Circuit
import pandas as pd
import matplotlib.pyplot as plt
from TORCphysics import visualization as vs

# TODO:
#  1.- Make the code work with stochastic topo binding.
#  2.- Run experiment to test that it works.
#  3.- Calibrate it so it works with Sam Meyers model, for now...

# Initial conditions
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'output'
frames = 100
series = True
continuation = False
tm = 'stochastic'
mm = 'uniform'
dt = 1.0
n_simulations = 1

for ns in range(n_simulations):

    # Load simulation
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt, tm, mm)

    my_circuit.name = my_circuit.name + '_' + str(ns)
    my_circuit.log.name = my_circuit.name
    my_circuit.print_general_information()
    my_circuit.run()

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 3

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'Raspberry': 'red'}
kwargs = {'linewidth': 2, 'ls': '--'}

fig, axs = plt.subplots(2, figsize=(width, 2 * height), tight_layout=True)

# Plot site response curves
# ---------------------------------------------------------
ax = axs[0]
to_ignore = [site.name for site in my_circuit.site_list if site.name != 'sum']
vs.plot_site_response_curves(my_circuit, ax)  #  , ignore=to_ignore)
# Let's ignore the

plt.savefig('calibration.png')