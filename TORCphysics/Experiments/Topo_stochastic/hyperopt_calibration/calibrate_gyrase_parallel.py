from TORCphysics.Experiments.Topo_stochastic.hyperopt_calibration.calibration_tools.calibration_tools \
    import run_single_stochastic
from TORCphysics import Circuit
from hyperopt import tpe, hp, fmin
import numpy as np
import concurrent.futures
from multiprocessing import Pool

# import sys

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# The idea is to calibrate topo_k_cat using the hyperopt library

# ----------------------------------------------------------------------------------------------------------------------
# INITIAL CONDITIONS
# ----------------------------------------------------------------------------------------------------------------------
# Optimization conditions
num_workers = 12
chunksize = 8
nsim = 96  # in total
ntests = 200  # number of tests for parametrization
# coutput = 'gyra_k_cat.txt'
coutput = 'gyra_kcat-alpha_kon.txt'
coutput = 'gyra_kcat-alpha_kon_koff.txt'
coutput = 'gyra_cat_kon_sigmoid.txt'
coutput = 'gyra_ks_sigmoid.txt'
k_cat_min = -20.0  # Ranges to vary k_cat
k_cat_max = -5.0
alpha_min = 0.5
alpha_max = 3.0
k_off_min = 0.1
k_off_max = 1.0
k_on_min = 0.002
k_on_max = 0.08
width_min = 0.005
width_max = 0.03
threshold_min = 0.006
threshold_max = 0.012

my_var = 'gyra_k_cat'
# my_vars = ['gyra_k_cat', 'alpha']
my_vars = ['k_cat', 'alpha', 'k_on', 'k_off']
my_vars = ['k_cat', 'k_on', 'width', 'threshold']
my_vars = ['k_cat', 'k_on', 'k_off', 'width', 'threshold']

# Some initial conditions
topo_k_on_0 = .005  # 0.0075
topo_k_cat_0 = -20.0
topo_concentration_0 = 0.0
gyra_k_on_0 = .005  # 0.0075
gyra_k_cat_0 = -20
gyra_concentration_0 = 0.25
gyra_k_off_0 = 0.5
topo_k_off_0 = 0.5

# Circuit conditions
circuit_filename_0 = 'circuit_3000bp_positive.csv'
sites_filename_0 = 'sites.csv'
enzymes_filename_0 = 'enzymes.csv'
environment_continuum_filename_0 = 'environment_continuum.csv'
environment_stochastic_filename_0 = 'environment_stochastic.csv'
tm = 'stochastic'
output_prefix = 'output'
frames = 2000
series = True
continuation = False
mm = 'uniform'
dt = .5


# ----------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------
def set_items(my_gyra_k_cat, my_gyra_k_on, my_gyra_k_off, my_gyra_width, my_gyra_threshold):
    items = []
    # for n in range(num_workers):
    for n in range(nsim):
        item = {
            'topo_concentration': topo_concentration_0,
            'topo_k_on': topo_k_on_0,
            'topo_k_cat': topo_k_cat_0,
            'topo_k_off': topo_k_off_0,
            'gyra_concentration': gyra_concentration_0,
#            'gyra_concentration': my_gyra_alpha * gyra_concentration_0,
            'gyra_k_on': my_gyra_k_on,
            'gyra_k_cat': my_gyra_k_cat,
            'gyra_width': my_gyra_width,
            'gyra_threshold': my_gyra_threshold,
            'gyra_k_off': my_gyra_k_off,
            'circuit_filename': circuit_filename_0,
            'sites_filename': sites_filename_0,
            'enzymes_filename': enzymes_filename_0,
            'environment_filename': environment_stochastic_filename_0,
            'output_prefix': output_prefix,
            'frames': frames,
            'series': series,
            'continuation': continuation,
            'dt': dt,
            'tm': 'stochastic',
            'mm': mm
        }
        items.append(item)
    return items


def objective_ProcessPoolExecutor(params):
    my_gyra_k_cat = params[my_vars[0]]
    my_gyra_k_on = params[my_vars[1]]
    my_gyra_k_off = params[my_vars[2]]
    my_gyra_width = params[my_vars[3]]
    my_gyra_threshold = params[my_vars[4]]
    items = set_items(my_gyra_k_cat, my_gyra_k_on, my_gyra_k_off, my_gyra_width, my_gyra_threshold)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit tasks to executor
        results = list(executor.map(run_single_stochastic, items, chunksize=chunksize))

        # Retrieve the superhelical densities
        # supercoiling = [future.result() for future in concurrent.futures.as_completed(futures)]

    my_supercoiling = np.zeros((frames + 1, nsim))
    for i, sigma in enumerate(results):
        my_supercoiling[:, i] = sigma
    # meean = np.mean(my_supercoiling, axis=1)
    # stdd  = np.std(my_supercoiling, axis=1)
    my_objective = np.sum(np.square(np.mean(my_supercoiling, axis=1) - sigma_continuum))
    return my_objective


# ----------------------------------------------------------------------------------------------------------------------
# PROCESS
# ----------------------------------------------------------------------------------------------------------------------

# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------
tm = 'continuum'
continuum_circuit = Circuit(circuit_filename_0, sites_filename_0, enzymes_filename_0, environment_continuum_filename_0,
                            output_prefix, frames, series, continuation, dt, tm, mm)
continuum_circuit.environmental_list[1].concentration = gyra_concentration_0
continuum_circuit.run()

# Get global supercoiling responses from the continuum case
mask = continuum_circuit.sites_df['type'] == 'circuit'  # This one contains global superhelical density
sigma_continuum = continuum_circuit.sites_df[mask]['superhelical'].to_numpy()

# Optimization case
# ----------------------------------------------------------------------------------------------------------------------
space = {
    #    my_var: hp.uniform(my_var, k_cat_min, k_cat_max)
    my_vars[0]: hp.uniform(my_vars[0], k_cat_min, k_cat_max),
    my_vars[1]: hp.uniform(my_vars[1], k_on_min, k_on_max),
    my_vars[2]: hp.uniform(my_vars[2], k_off_min, k_off_max),
    my_vars[3]: hp.uniform(my_vars[3], width_min, width_max),
    my_vars[4]: hp.uniform(my_vars[4], threshold_min, threshold_max)
    #    my_vars[1]: hp.uniform(my_vars[1], alpha_min, alpha_max),
#    my_vars[3]: hp.uniform(my_vars[3], k_off_min, k_off_max)
}

# We can distribute tuning across our Spark cluster
# by calling `fmin` with a `SparkTrials` instance.
# sparktrials = SparkTrials(parallelism=4)

best = fmin(
    #  fn=objective_gyra_k_cat, # Objective Function to optimize
    #    fn=objective_gyra_k_cat_alpha_k_on_parallel_concurrent,  # Objective Function to optimize
    fn=objective_ProcessPoolExecutor,  # Objective Function to optimize
    space=space,  # Hyperparameter's Search Space
    algo=tpe.suggest,  # Optimization algorithm (representative TPE)
    max_evals=ntests  # Number of optimization attempts
)
print(best)
print(best.values())
best_result = np.zeros(5)
#  best_result[0] = best[my_var]  # Extract the best parameter
best_result[0] = best[my_vars[0]]  # Extract the best parameter
best_result[1] = best[my_vars[1]]
best_result[2] = best[my_vars[2]]
best_result[3] = best[my_vars[3]]
best_result[4] = best[my_vars[4]]
np.savetxt(coutput, best_result)


# ----------------------------------------------------------------------------------------------------------------------
# Previous functions that you should ignore...
def objective_gyra_k_cat(params):
    my_gyra_k_cat = params[my_var]
    supercoiling = run_many_stochastic(topo_concentration_0, topo_k_on_0, topo_k_cat_0,
                                       gyra_concentration_0, gyra_k_on_0, my_gyra_k_cat,
                                       circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                       environment_stochastic_filename_0, output_prefix, frames, series, continuation,
                                       dt, 'stochastic', mm)
    my_objective = np.sum(np.square(np.mean(supercoiling, axis=1) - sigma_continuum))
    return my_objective


def objective_gyra_k_cat_alpha_(params):
    my_gyra_k_cat = params[my_vars[0]]
    my_gyra_alpha = params[my_vars[1]]
    supercoiling = run_many_stochastic(topo_concentration_0, topo_k_on_0, topo_k_cat_0,
                                       my_gyra_alpha * gyra_concentration_0, gyra_k_on_0, my_gyra_k_cat,
                                       circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                       environment_stochastic_filename_0, output_prefix, frames, series, continuation,
                                       dt, 'stochastic', mm)
    my_objective = np.sum(np.square(np.mean(supercoiling, axis=1) - sigma_continuum))
    return my_objective


def objective_gyra_k_cat_alpha_k_on_parallel_multiprocessing(params):
    my_gyra_k_cat = params[my_vars[0]]
    my_gyra_alpha = params[my_vars[1]]
    my_gyra_k_on = params[my_vars[2]]
    items = set_items(my_gyra_k_cat, my_gyra_k_on, my_gyra_alpha)
    with Pool(processes=num_workers) as pool:
        supercoiling = pool.map(run_many_stochastic, items)

    my_supercoiling = np.zeros((frames + 1, nsim * num_workers))
    s = 0
    for i, sigma in enumerate(supercoiling):
        for n in range(nsim):
            my_supercoiling[:, s] = sigma[:, n]
            s += 1
    # meean = np.mean(my_supercoiling, axis=1)
    # stdd  = np.std(my_supercoiling, axis=1)
    my_objective = np.sum(np.square(np.mean(my_supercoiling, axis=1) - sigma_continuum))
    return my_objective


def objective_gyra_k_cat_alpha_k_on_parallel_concurrent(params):
    my_gyra_k_cat = params[my_vars[0]]
    my_gyra_alpha = params[my_vars[1]]
    my_gyra_k_on = params[my_vars[2]]
    items = set_items(my_gyra_k_cat, my_gyra_k_on, my_gyra_alpha)
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:

        # Submit tasks to executor
        futures = [executor.submit(run_many_stochastic, item) for item in items]

        # Retrieve the superhelical densities
        supercoiling = [future.result() for future in concurrent.futures.as_completed(futures)]

    my_supercoiling = np.zeros((frames + 1, nsim * num_workers))
    s = 0
    for i, sigma in enumerate(supercoiling):
        for n in range(nsim):
            my_supercoiling[:, s] = sigma[:, n]
            s += 1
    # meean = np.mean(my_supercoiling, axis=1)
    # stdd  = np.std(my_supercoiling, axis=1)
    my_objective = np.sum(np.square(np.mean(my_supercoiling, axis=1) - sigma_continuum))
    return my_objective
