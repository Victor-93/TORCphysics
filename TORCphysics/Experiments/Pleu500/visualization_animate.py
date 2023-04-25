from TORCphysics import Circuit
import pandas as pd
from TORCphysics import visualization as vs


# Description
# ---------------------------------------------------------
# This script is produced with the intention to test the visualization module.
# Even though this is not a formal test, we use it to check that we are producing the desired animations.

# Inputs
# ---------------------------------------------------------
csites_df = 'Pleu500_0_sites_df.csv'
cenzymes_df = 'Pleu500_0_enzymes_df.csv'
cenvironment_df = 'Pleu500_0_environment_df.csv'
sites_df = pd.read_csv(csites_df, sep=',')
enzymes_df = pd.read_csv(cenzymes_df, sep=',')

log_file = 'Pleu500_0.log'

colors_dict = {'tetA': 'yellow', 'CDS': 'green', 'mKalama1': 'blue', 'Raspberry': 'red'}

circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'output'
frames = 1000
series = True
continuation = False
tm = 'continuum'
mm = 'uniform'
dt = 1.0
n_simulations = 1
bridge_time = 40000

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt, tm, mm)

output = 'test_animation'
out_format = '.gif'
vs.create_animation_linear(my_circuit, sites_df, enzymes_df, output, out_format,
                           site_type='gene', site_colours=colors_dict)
