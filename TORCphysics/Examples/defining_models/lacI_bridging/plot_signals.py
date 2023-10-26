from TORCphysics import Circuit
from TORCphysics.src import visualization as vs
import pandas as pd
import matplotlib.pyplot as plt

# Inputs
# ---------------------------------------------------------
csites_df = 'test1__sites_df.csv'
cenzymes_df = 'test1__enzymes_df.csv'
cenvironment_df = 'test1__environment_df.csv'
sites_df = pd.read_csv(csites_df, sep=',')
enzymes_df = pd.read_csv(cenzymes_df, sep=',')

log_file = 'test1__.log'

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'Raspberry': 'red',
               'lac1': 'white', 'lac2': 'white'}

circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'out'
frames = 2000
series = True
continuation = False
dt = .5

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)
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
dt = 0.5  # This should be extracted from the log file
# Create Figure
# ---------------------------------------------------------
fig, axs = plt.subplots(1, figsize=(width,  height), tight_layout=True)

# Sites rate curves - Let's plot the rates modulated by supercoiling
# ---------------------------------------------------------
ax = axs
vs.plot_signal_profiles(my_circuit, sites_df, ax, site_type='gene', ignore='CDS', colors=colors_dict)
#vs.plot_signal_profiles(my_circuit, sites_df, ax, ignore=['ori', 'test1'], colors=colors_dict)
ax.set_title('Signals')
plt.savefig('signals.png')

