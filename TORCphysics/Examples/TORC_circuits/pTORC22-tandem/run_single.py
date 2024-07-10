from TORCphysics import Circuit

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
dt = 1.

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)
my_circuit.print_general_information()
my_circuit.run()
