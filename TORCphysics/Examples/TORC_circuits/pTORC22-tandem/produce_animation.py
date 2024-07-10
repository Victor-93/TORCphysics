from TORCphysics import Circuit
import pandas as pd
from TORCphysics import visualization as vs

# Description
# ---------------------------------------------------------
# Here, we will be testing the visualization package. We want it to have it ready for a conference

# TODO: 1.- Try using objects for topos (clouds) - this seems more complicated, let's leave it for later.
#  2.- Use buckets that fills with transcripts.
#  3.- Make nice visualization (maybe use params from k_on fixed? of topos).
#  4.- Maybe you could make a linear version of the plasmid, right?
# Inputs
# ---------------------------------------------------------
with_topoI = True

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = '../enzymes.csv'
if with_topoI:
    environment_filename = '../environment.csv'
    output_prefix = 'WT_system'
else:
    environment_filename = '../environment_notopoI.csv'
    output_prefix = 'notopoI_system'
frames = 250 #50000
series = True
continuation = False
dt = 1.0

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
colors_dict = {'tetA': '#d8b41a', 'rop': 'silver', 'mKalama1': '#0051ff', 'mRaspberry': '#e30000', 'antitet': 'orange',
               'bla':'purple'}

output = 'animation_ON_30fps'
out_format = '.mp4'

#vs.create_animation_linear(my_circuit, sites_df, enzymes_df, output, out_format,
#                           site_type='gene', site_colours=colors_dict)
vs.create_animation_linear_artist(my_circuit, sites_df, enzymes_df, output, out_format,
                           site_type='gene', site_colours=colors_dict)
