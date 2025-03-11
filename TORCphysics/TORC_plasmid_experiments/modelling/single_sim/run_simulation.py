from TORCphysics import Circuit, params
import pandas as pd
from TORCphysics import visualization as vs
import numpy as np

# Description
# ----------------------------------------------------------------------------------------------------------------------
# This script simply runs TORCPhysiscs for the TORC Plasmid

# Simulation conditions
# ----------------------------------------------------------------------------------------------------------------------
dt = 1.0
initial_time = 0
final_time = 5000
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Circuit initial conditions
# ----------------------------------------------------------------------------------------------------------------------
enzymes_filename = None  # Because nothing is bound initially
series = True
continuation = False

# Minimum system conditions - Choose this for the minimal promoter system
# The disp-sites* are the ones for which lacI represses PleuWT initiation
circuit_filename = '../circuit_min-linear.csv'
sites_filename = 'disp-sites_min-linear.csv'
#output_prefix = 'minimum_Ecoli'
output_prefix = 'minimum_Sal'
#output_prefix = 'no_topos'

# Complete system conditions - Choose for the complete promoter
#circuit_filename = '../circuit_complete-linear.csv'
#sites_filename = 'disp-sites_min-linear.csv'
#output_prefix = 'complete_Ecoli'

# Environmentals
#environment_filename = 'Ecoli_environment.csv' # For Ecoli
environment_filename = 'Salmonella_environment.csv' # For Salmonella
#environment_filename = 'no_topos_environment.csv'

# Simulation
# ----------------------------------------------------------------------------------------------------------------------
# Load base version of circuit with all sites, environmentals, but without models
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)
my_circuit.print_general_information()

#enzymes_df, sites_df, environmental_df = my_circuit.run_return_dfs()

my_circuit.run()
