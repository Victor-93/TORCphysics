from TORCphysics import Circuit, Enzyme
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import pandas as pd

# TODO:
#  1.- Make the code work with stochastic topo binding.
#  2.- Run experiment to test that it works.
#  3.- Calibrate it so it works with Sam Meyers model, for now...

# Initial conditions
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'output'
frames = 1000
series = True
continuation = False
tm = 'stochastic'
mm = 'uniform'
dt = 1.0
n_simulations = 1

for ns in range(n_simulations):
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt, tm, mm)

    my_circuit.name = my_circuit.name + '_' + str(ns)
    my_circuit.log.name = my_circuit.name
    my_circuit.print_general_information()
    my_circuit.run()
