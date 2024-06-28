import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from TORCphysics import Circuit
from TORCphysics import visualization as vs


# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by comp_feed.py.
names = ['weak', 'medium', 'strong']

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = 'test_circuit.csv'
sites_filename = 'test_sites.csv'
enzymes_filename = None
environment_filename = '../torcphys/environment.csv'
frames = 5000 #50000
series = True
continuation = False
dt = 0.25

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     'output_prefix', frames, series, continuation, dt)

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4

colors_dict = {'weak': 'green', 'medium': 'blue', 'strong': 'red'}
kwargs = {'linewidth': 2, 'ls': '-'}

# Plot responses
# -------------------------------------
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
fig.suptitle(' rates')
vs.plot_site_response_curves(my_circuit, axs=axs, colors=colors_dict, site_type='gene')
plt.savefig('responses.png')
plt.show()