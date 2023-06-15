import concurrent.futures
from TORCphysics import parallelization_tools as pt
# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# We want to run a multiple simulations of the Pleu500 with stochastic topoisomerase activities

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------

# Paralelization conditions
n_workers = 12
n_simulations = 96  # in total
simulations_per_worker = n_simulations // n_workers

# Circuit conditions
circuit_filename = '../../../circuit.csv'
sites_filename = 'sites_maxmin.csv'
enzymes_filename = '../../../enzymes_OFF.csv'
environment_filename = 'environment_stochastic.csv'
output_prefix = 'out'
frames = 8000
series = True
continuation = False
tm = 'stochastic'
mm = 'uniform'
dt = .5

# Create a ProcessPoolExecutor
with concurrent.futures.ProcessPoolExecutor() as executor:
    # Submit simulation tasks to the executor
    futures = []
    for worker_id in range(n_workers):
        start_simulation = worker_id * simulations_per_worker + 1
        end_simulation = (worker_id + 1) * simulations_per_worker
        for simulation_number in range(start_simulation, end_simulation + 1):
            item = pt.set_item(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                               output_prefix, frames, series, continuation, dt, tm, mm, simulation_number)
            future = executor.submit(pt.run_single_simulation, item)
            futures.append(future)
#            print(worker_id,simulation_number)

    # Wait for all simulations to complete
    concurrent.futures.wait(futures)
