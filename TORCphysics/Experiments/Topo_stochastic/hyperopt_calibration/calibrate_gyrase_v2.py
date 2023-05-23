from TORCphysics.Experiments.Topo_stochastic.hyperopt_calibration.calibration_tools.calibration_tools \
    import run_single_stochastic
from TORCphysics import Circuit
from hyperopt import tpe, hp, fmin
import numpy as np
import concurrent.futures

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
ntests = 100  # number of tests for parametrization
coutput = 'gyra_calibration.txt'
k_cat_min = -20.0  # Ranges to vary k_cat
k_cat_max = -5.0
alpha_min = 0.1
alpha_max = 3.0
width_min = 0.001
width_max = 0.03
threshold_min = 0.001
threshold_max = 0.02

k_off = 0.5
k_on = 0.005
concentration = 0.25

my_vars = ['k_cat', 'alpha', 'width', 'threshold']

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
def set_items(topo_kcat, topo_kon, topo_koff, topo_alpha, topo_concentration, topo_width, topo_threshold,
              gyra_kcat, gyra_kon, gyra_koff, gyra_alpha, gyra_concentration, gyra_width, gyra_threshold):
    items = []
    # for n in range(num_workers):
    for n in range(nsim):
        item = {
            'topo_concentration': topo_concentration,
            'topo_k_on': topo_kon,
            'topo_k_cat': topo_kcat,
            'topo_k_off': topo_koff,
            'topo_alpha': topo_alpha,
            'topo_width': topo_width,
            'topo_threshold': topo_threshold,
            'gyra_concentration': gyra_concentration,
            'gyra_k_on': gyra_kon,
            'gyra_k_off': gyra_koff,
            'gyra_k_cat': gyra_kcat,
            'gyra_alpha': gyra_alpha,
            'gyra_width': gyra_width,
            'gyra_threshold': gyra_threshold,
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
    my_kcat = params[my_vars[0]]
    my_alpha = params[my_vars[1]]
    my_width = params[my_vars[2]]
    my_threshold = params[my_vars[3]]
    # For Gyrase
    items = set_items( topo_kcat=0.0, topo_kon=k_on, topo_koff=k_off, topo_alpha=0.0,
                       topo_concentration=0.0, topo_width=0.1, topo_threshold=0.0,
                       gyra_kcat=my_kcat, gyra_kon=k_on, gyra_koff=k_off, gyra_alpha=my_alpha,
                       gyra_concentration=concentration, gyra_width=my_width, gyra_threshold=my_threshold)
    # For Topo
    #  items = set_items( topo_kcat=my_kcat, topo_kon=k_on, topo_koff=k_off, topo_alpha=my_alpha,
    #                   topo_concentration=concentration, topo_width=my_width, topo_threshold=my_threshold,
    #                   gyra_kcat=0.0, gyra_kon=k_on, gyra_koff=k_off, gyra_alpha=0.0,
    #                   gyra_concentration=0.0, gyra_width=0.1, gyra_threshold=0.0)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit tasks to executor
        results = list(executor.map(run_single_stochastic, items, chunksize=chunksize))

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
    my_vars[0]: hp.uniform(my_vars[0], k_cat_min, k_cat_max),
    my_vars[1]: hp.uniform(my_vars[1], alpha_min, alpha_max),
    my_vars[2]: hp.uniform(my_vars[2], width_min, width_max),
    my_vars[3]: hp.uniform(my_vars[3], threshold_min, threshold_max)
}

best = fmin(
    fn=objective_ProcessPoolExecutor,  # Objective Function to optimize
    space=space,  # Hyperparameter's Search Space
    algo=tpe.suggest,  # Optimization algorithm (representative TPE)
    max_evals=ntests  # Number of optimization attempts
)
print(best)
print(best.values())
best_result = np.zeros(4)
best_result[0] = best[my_vars[0]]  # Extract the best parameter
best_result[1] = best[my_vars[1]]
best_result[2] = best[my_vars[2]]
best_result[3] = best[my_vars[3]]
np.savetxt(coutput, best_result)
