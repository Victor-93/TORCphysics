import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from TORCphysics import Circuit
from TORCphysics import visualization as vs

# TODO ALGO FALTA PARA LAS RESPONSES

# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by comp_feed.py.
names = ['weak', 'medium', 'strong']

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = 'test_circuit.csv'
sites_filename = 'test_sites.csv'
enzymes_filename = None
environment_filename = 'environment.csv'
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

yfp_c = 'gold'
bla_c = 'purple'
tetA_c = 'green'
antitet_c = 'cyan'

colors_dict = {'antitet': 'cyan', 'bla': 'purple', 'tetA': 'green', 'PleuWT': 'gold'}
kwargs = {'linewidth': 2, 'ls': '-'}

# Plot responses
# -------------------------------------
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
fig.suptitle(' rates')
vs.plot_site_response_curves(my_circuit, axs=axs, colors=colors_dict, site_type='gene')
plt.savefig('responses.png')
plt.show()