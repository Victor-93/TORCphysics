import sys
#sys.path.append("/users/ph1vvb")
import numpy as np
import pandas as pd
from hyperopt import tpe, hp, fmin, Trials
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import pickle



# Description
# --------------------------------------------------------------
# Following same logic than calibrate_tracking_nsets, this script runs the calibration using the function:
# single_case_RNAPTracking_calibration_nsets_p2()
# This function uses a different parallelization strategy to speed up the calibration process.
# It does this by doing two parallelization pools.

# This script launches the calibration process for the TopoIRNAPTracking model.
# It uses the parameters and models of the environment.csv file within this folder
# Runs multiple sets of simulations, to obtain KDEs and then average them to smooth the curve and obtain more
# stable correlations.

# Parallelization conditions
# --------------------------------------------------------------
# NOTE: The number of sets or KDEs that are obtain are n_subsets * n_sets.
#       The number of simulations launched per n_subset is n_inner_workers
#       In total, the number of simulations launched is n_subset * n_sets * n_inner_workers
# WARNING: It is very important that the number of n_workers is within the capabilities of your system.
#          So choose carefully n_workers and n_sets
n_workers = 64#12  # Total number of workers (cpus)
n_sets = 7#3  # Number of outer sets
n_subsets = 4#7#2  # Number of simulations per set
n_inner_workers = n_workers // (n_sets+1)  # Number of workers per inner pool
                                           # +1 because one worker is spent in creating the outer pool

tests = 200  #60  # number of tests for parametrization

print('Doing parallelization process for:')
print('n_workers', n_workers)
print('n_sets', n_sets)
print('n_subsets', n_subsets)
print('n_inner_workers', n_inner_workers)
print('hyperopt tests', tests)
print('Total number of simulations per test:', n_sets * n_subsets * n_inner_workers)
print('Total number of actual workers:', n_sets * (1 + n_inner_workers))

# Simulation conditions
# --------------------------------------------------------------
file_out = 'calibration_RNAPTracking_nsets_p2'
dt = 0.25
initial_time = 0
final_time = 1000#500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Reference - It is the reference density of topos when there is no gene that we will use to calculate the
#             fold enrichment.
# --------------------------------------------------------------
reference_path = '../noRNAP/'  # Path with the reference kdes
reference_RNAP = '../RNAP_TU.txt'

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = '../sites.csv'
enzymes_filename = None  #'../enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'topoIRNAPtrack-uniform'
series = True
continuation = False

# Let's asume it starts at a value where topos are equilibrated, so we assume the steady state.

# Models to calibrate
# -----------------------------------
# Topoisomerase I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'TopoIRecognitionRNAPTracking'
topoI_params_csv = '../calibration_topoI.csv'  # These are parameters for the Recognition curve

# RANGES FOR RANDOM SEARCH
# -----------------------------------
RNAP_dist_min = 20
RNAP_dist_max = 500
fold_change_min = .1
fold_change_max = 50

# TARGETS FOR OPTIMIZATION
# -----------------------------------
target_FE = 1.38 #1.68  # Target fold-enrichment
target_CO = 1.0  # Target correlation between topo I and RNAP densities.
x_spacing = 10.0  # The spacing I want at the moment of doing interpolation.

# nbins is the number of bins to use when calculating the kde
target_dict = {'target_FE': target_FE, 'target_CO': target_CO, 'target_gene': 'reporter',
               'enzymes_names': ['RNAP', 'topoI', 'gyrase']}


# ----------------------------------------------------------------------------------------------------------------------
# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------

# If calibrating=True, then only the objective function is returned (used by hyperopt). If it is false, then it returns
# the outputs
def objective_function(params, calibrating=True):
    # We need to prepare the inputs.
    # At the moment, we only have one system.

    # Global dictionaries
    # ------------------------------------------
    global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                   'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                   'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                   'frames': frames, 'dt': dt,
                   'n_workers': n_workers, 'n_sets': n_sets, 'n_subsets': n_subsets, 'n_inner_workers': n_inner_workers,
                   'DNA_concentration': 0.0}

    # Variation dictionaries
    # ------------------------------------------

    # Topoisomerase I
    name = topoI_name
    object_type = topoI_type
    binding_model_name = topoI_binding_model_name
    binding_oparams = {'k_on': float(topoI_params['k_on'][0]), 'width': float(topoI_params['width'][0]),
                       'threshold': float(topoI_params['threshold'][0]),
                       'RNAP_dist': params['RNAP_dist'],
                       'fold_change': params['fold_change']}

    topoI_variation = {'name': name, 'object_type': object_type,
                       'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams}

    # Create lists of conditions for each system
    # ------------------------------------------

    # Global dictionaries
    global_dict_list = [global_dict]

    # List of lists of variations
    variations_list = [[topoI_variation]]

    # Arrays with position densities to calculate fold change
    list_reference = [reference_dict]

    # Finally, run objective function.
    # ------------------------------------------
    my_objective, output_dict = tct.single_case_RNAPTracking_calibration_nsets_2scheme(global_dict_list, variations_list,
                                                                               list_reference,
                                                                               target_dict)

#    return my_objective
    if calibrating:
        return my_objective
    else:
        return my_objective, output_dict


# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------

# Let's load the circuit, so we can extract some information
# -----------------------------------
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

# Get target site
target_gene = [site for site in my_circuit.site_list if site.name == target_dict['target_gene']][0]
RNAP_env = [environment for environment in my_circuit.environmental_list if environment.name == 'RNAP'][0]

# Define x-axes
x_system = tct.get_interpolated_x(1, my_circuit.size)
x_gene = tct.get_interpolated_x(target_gene.start - RNAP_env.size, target_gene.end)

reference_dict = {}
for name in target_dict['enzymes_names']:  # Reference topos

    # Load reference file
    if name != 'RNAP':
        kde_ref = np.loadtxt(reference_path + 'reference_' + name + '.txt')

        # And get interpolated data - I do it here, so we don't have to do it again in the parallelization and save some
        # time. It is saved to the reference_dict by the way
        reference_dict[name] = tct.get_interpolated_kde(kde_ref[:, 0], kde_ref[:, 1], x_system)

#Reference RNAP
RNAP_TU = np.loadtxt(reference_RNAP)
length_TU = len(RNAP_TU)
dx = (target_gene.end - (target_gene.start - RNAP_env.size))/length_TU
TU_x = np.arange(target_gene.start - RNAP_env.size, target_gene.end, dx)
reference_dict['RNAP'] = tct.get_interpolated_kde(TU_x, RNAP_TU, x_gene)

# Load topoI params
topoI_params = pd.read_csv(topoI_params_csv).to_dict()

# Optimization
# -----------------------------------------------------
trials = Trials()
space = {
    # Topo I params
    'RNAP_dist': hp.uniform('RNAP_dist', RNAP_dist_min, RNAP_dist_max),
    'fold_change': hp.uniform('fold_change', fold_change_min, fold_change_max)
}

# Save the current standard output
original_stdout = sys.stdout
# Define the file where you want to save the output
output_file_path = file_out + '.info'

# Open the file in write mode
with open(output_file_path, 'w') as f:
    # Redirect the standard output to the file
    sys.stdout = f  ##

    # Your code that prints to the screen
    print("Hello, this is the info file for the calibration of Topo I Tracking RNAP.")
    print("Topo I Binding Model = " + topoI_binding_model_name)
    print('n_workers', n_workers)
    print('n_sets', n_sets)
    print('n_subsets', n_subsets)
    print('n_inner_workers', n_inner_workers)
    print('hyperopt tests', tests)
    print('Total number of simulations per test:', n_sets * n_subsets * n_inner_workers)
    print('Total number of actual workers:', n_sets * (1 + n_inner_workers))

    best = fmin(
        fn=objective_function,  # Objective Function to optimize
        space=space,  # Hyperparameter's Search Space
        algo=tpe.suggest,  # Optimization algorithm (representative TPE)
        max_evals=tests,  # Number of optimization attempts
        trials=trials
    )

    print(" ")
    print("Optimal parameters found from random search: ")
    print(best)

best_df = pd.DataFrame.from_dict([best])

# Save model!
# --------------------------------------------------------------------------
k_on = float(topoI_params['k_on'][0])
width = float(topoI_params['width'][0])
threshold = float(topoI_params['threshold'][0])

# This one will have all the values of the RNAPTracking Model
complete_df = pd.DataFrame(columns=['k_on', 'width', 'threshold', 'RNAP_dist', 'fold_change'])
complete_df['RNAP_dist'] = best_df['RNAP_dist']
complete_df['fold_change'] = best_df['fold_change']
complete_df['k_on'] = k_on
complete_df['width'] = width
complete_df['threshold'] = threshold
complete_df.to_csv(file_out + '.csv', index=False, sep=',')

# Let's save trials info (params and loses)
# --------------------------------------------------------------------------
params_df = pd.DataFrame(columns=['test', 'loss', 'RNAP_dist', 'fold_change'])
for n in range(tests):
    tdi = trials.trials[n]  # dictionary with results for test n
    lo = trials.trials[n]['result']['loss']  # loss
    va = trials.trials[n]['misc']['vals']  #values
    # Add a new row using append method
    new_row = pd.DataFrame({
        'test': n, 'loss': lo, 'RNAP_dist': va['RNAP_dist'], 'fold_change': va['fold_change'],
        'k_on': k_on, 'width': width, 'threshold': threshold
    })
    #    params_df.append(new_row, ignore_index=True)
    params_df = pd.concat([params_df, new_row], ignore_index=True)

params_df.to_csv(file_out+'-values.csv', index=False, sep=',')

# Let's run the function once more with the best params to produce the data so we can then just plot it.
# --------------------------------------------------------------------------
objective, output_dict = objective_function(params=best, calibrating=False)

# Save the dictionary to a file
with open(file_out+'.pkl', 'wb') as file:
    pickle.dump(output_dict, file)