import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from TORCphysics import Circuit
from TORCphysics import visualization as vs


# Description
# ---------------------------------------------------------
# I will process and analyse the simulations produced by comp_feed.py.
names = ['circuit', 'tetA', "mKalama1", 'mRaspberry', 'antitet', 'bla']
n_simulations = 12
with_topoI = False

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = None
if with_topoI:
    environment_filename = '../environment.csv'
    output_prefix = 'WT_system'
    enzyme_names = ['RNAP', 'topoI', 'gyrase']
else:
    environment_filename = '../environment_notopoI.csv'
    output_prefix = 'notopoI_system'
    enzyme_names = ['RNAP', 'gyrase']
frames = 5000 #50000
series = True
continuation = False
dt = 0.25

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

# Figure initial conditions
# ---------------------------------------------------------
width = 8
height = 4

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'mRaspberry': 'red', 'lac1': 'green',
               'lac2': 'green', 'antitet':'green', 'circuit':'black', 'RNAP': 'black', 'topoI': 'red', 'gyrase':'cyan',
               'bla': 'purple'}
kwargs = {'linewidth': 2, 'ls': '-'}


names_genes = ['tetA', 'mKalama1', 'antitet', 'mRaspberry', 'bla']


# Plot responses
# -------------------------------------
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
fig.suptitle(my_circuit.name + ' rates')
vs.plot_site_response_curves(my_circuit, axs=axs, colors=colors_dict, site_type='gene')
plt.savefig('responses.png')