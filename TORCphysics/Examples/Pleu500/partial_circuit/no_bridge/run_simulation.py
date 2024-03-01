from TORCphysics import Circuit


# Let's initialize circuit
circuit_filename = '../circuit.csv'
sites_filename = '../sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'nobridge'
frames = 10000
series = True
continuation = False
dt = 1.0#0.25
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)
my_circuit.print_general_information()
my_circuit.run()

