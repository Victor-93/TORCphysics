from TORCphysics import Circuit
import pandas as pd
from TORCphysics import visualization as vs

# Description
# ---------------------------------------------------------
# I just want to visualize the simulation with stochastic topos and stochastic transcription

# Inputs
# ---------------------------------------------------------
csites_df = 'Pleu500_visual_sites_df.csv'
cenzymes_df = 'Pleu500_visual_enzymes_df.csv'
cenvironment_df = 'Pleu500_visual_environment_df.csv'
sites_df = pd.read_csv(csites_df, sep=',')
enzymes_df = pd.read_csv(cenzymes_df, sep=',')

log_file = 'Pleu500_out.log'

colors_dict = {'tetA': '#d8b41a', 'CDS': 'silver', 'mKalama1': '#0051ff', 'Raspberry': '#e30000'}

circuit_filename = '../../../circuit.csv'
sites_filename = 'sites_maxmin.csv'
enzymes_filename = '../../../enzymes.csv'
environment_filename = 'environment_stochastic.csv'
output_prefix = 'out'
frames = 1000
series = True
continuation = False
tm = 'stochastic'
mm = 'uniform'
dt = .25

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt, tm, mm)

output = 'animation'
out_format = '.gif'

vs.create_animation_linear(my_circuit, sites_df, enzymes_df, output, out_format,
                           site_type='gene', site_colours=colors_dict)
