from TORCphysics import Circuit

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# We want to run a single simulation of the Pleu500 with stochastic topoisomerase activities

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
circuit_filename = '../../../circuit.csv'
sites_filename = 'sites_maxmin.csv'
#sites_filename = 'sites_sam.csv'
enzymes_filename = '../../../enzymes_OFF.csv'
environment_filename = 'environment_stochastic.csv'
output_prefix = 'out'
frames = 8000
series = True
continuation = False
tm = 'stochastic'
mm = 'uniform'
dt = .5

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt, tm, mm)
my_circuit.print_general_information()
my_circuit.run()

