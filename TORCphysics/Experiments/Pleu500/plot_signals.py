from TORCphysics import analysis as an
import pandas as pd
import matplotlib.pyplot as plt

# Description
# ---------------------------------------------------------
# This is script is intended as an example to show how to process the outputs of TORCphysics.
# We will plot: 1.- the binding/transcription signals. 2.- Cross-correlations. 3.- Binding curves.
# 4.- Global supercoiling & supercoiling at site. 5.- topoisomerase activity curves (continuum). 6.- Rates.
# 7.- Supercoiling distributions at sites.

# Inputs
# ---------------------------------------------------------
csites_df = 'Pleu500_0_sites_df.csv'
cenzymes_df = 'Pleu500_0_enzymes_df.csv'
cenvironment_df = 'Pleu500_0_environment_df.csv'

sites_csv = 'Pleu500_0_sites_output.csv'
environment_csv = 'Pleu500_0_environment_output.csv'

log_file = 'Pleu500_0.log'

# Figure initial conditions
# ---------------------------------------------------------
width = 6
height = 3

# Better use these for the colors of genes...
# Sort them according the input file...
colors = []
colors.append("yellow")
colors.append("black")
colors.append("blue")
colors.append("magenta")
colors.append("red")
colors.append("green")
colors.append("cyan")
colors.append("black")

# Load inputs
# ---------------------------------------------------------
sites_df = pd.read_csv( csites_df, sep=',')
dt = 1.0  # This should be extracted from the log file
# Create Figure
# ---------------------------------------------------------
fig, axs = plt.subplots(2, figsize=(width, 2*height), tight_layout=True)

# Signals
# ---------------------------------------------------------
ax = axs[0]

signals, names = an.build_signal_by_type(sites_df, 'gene')
time = signals()
for i, signal in enumerate(signals):
    ax.plot()
