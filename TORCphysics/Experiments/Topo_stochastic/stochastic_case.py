from TORCphysics import Circuit
from TORCphysics import visualization as vs

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# The purpose of this script is to spead up a little bit the code, by looking at what features make it slower
# ----------------------------------------------------------------------------------------------------------------------

# Circuit conditions
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
frames = 1000
series = True
continuation = False
mm = 'uniform'
dt = 1.0

# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------
tm = 'stochastic'
output_prefix = 'stochastic'
environment_filename = 'environment.csv'
n_simulations = 1

# Load simulation
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt, tm, mm)
my_circuit.print_general_information()
my_circuit.run()
