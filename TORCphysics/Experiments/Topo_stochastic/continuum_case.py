from TORCphysics import Circuit

# Initial conditions
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment_continuum.csv'
output_prefix = 'continuum'
frames = 1000
series = True
continuation = False
tm = 'continuum'
mm = 'uniform'
dt = 1.0
n_simulations = 1

# Load simulation
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt, tm, mm)

my_circuit.print_general_information()
my_circuit.run()
