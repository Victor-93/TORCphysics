import multiprocessing
from TORCphysics import parallelization_tools as pt
import sys


# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# Uses multiprocessing instead of concurrent.futures
# We want to run a multiple simulations of the Pleu500 with stochastic topoisomerase activities

def run_single_simulation(simulation_number):
    item = pt.set_item(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                       output_prefix, frames, series, continuation, dt, tm, mm, simulation_number)
    pt.run_single_simulation(item)


# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------

# Paralelization conditions
#n_workers = 2  # 12
n_simulations = 100  # 96  # in total
#simulations_per_worker = n_simulations // n_workers

# Circuit conditions
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

# Create a multiprocessing pool
pool = multiprocessing.Pool()
items = []
for n in range(n_simulations):
    item = pt.set_item(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                       output_prefix, frames, series, continuation, dt, tm, mm, n)
    items.append(item)
pool.map(pt.run_single_simulation, items)



