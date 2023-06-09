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
fig, axs = plt.subplots(1, figsize=(width, 2 * height), tight_layout=True)

# Sites rate curves - Let's plot the rates modulated by supercoiling
# ---------------------------------------------------------
ax = axs
i = -1
vs.plot_site_response_curves(my_circuit, ax)
ax.set_title('Response curves')
#for site in my_circuit.site_list:
#    if site.site_type == 'gene':
#        i += 1
#        rate, x = an.site_activity_curves(site)
#        ax.plot(x, rate, color=colors[i], label=site.name)

#ax_params(axis=ax, xl=r'\sigma', yl=r'Initiation rate ($s^{-1}$)', grid=True, legend=True)
plt.show()


plt.savefig('responses.png')
