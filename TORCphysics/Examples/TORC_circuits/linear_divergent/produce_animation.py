from TORCphysics import Circuit
import pandas as pd
from TORCphysics import visualization as vs

# Description
# ---------------------------------------------------------
# Here, we will be testing the visualization package. We want it to have it ready for a conference

# Inputs
# ---------------------------------------------------------
with_topoI = True

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = None
if with_topoI:
    environment_filename = '../environment.csv'
    output_prefix = 'WT_system'
else:
    environment_filename = '../environment_notopoI.csv'
    output_prefix = 'notopoI_system'
frames = 2000#500 #50000
series = True
continuation = False
dt = 0.5

# Load circuit and dataframes
# --------------------------------------------------------------
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)


csites_df = my_circuit.name + '_' + output_prefix + '_sites_df.csv'
cenzymes_df = my_circuit.name + '_' + output_prefix + '_enzymes_df.csv'
cenvironment_df = my_circuit.name + '_' + output_prefix + '_environment_df.csv'
sites_df = pd.read_csv(csites_df, sep=',')
enzymes_df = pd.read_csv(cenzymes_df, sep=',')

# Animation stuff
# --------------------------------------------------------------
colors_dict = {'tetA': '#d8b41a', 'mKalama1': '#0051ff', 'mRaspberry': '#e30000',
               'bla': '#ca07e8'}

output = 'animation_ON_30fps'
out_format = '.mp4'

vs.create_animation_linear_artist(my_circuit, sites_df, enzymes_df, output, out_format,
                           site_type='gene', site_colours=colors_dict)
