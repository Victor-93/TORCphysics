import numpy as np
import pandas as pd
from hyperopt import tpe, hp, fmin, Trials
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
import sys

# Description
# --------------------------------------------------------------
# This script launches the calibration process for the TopoIRNAPTracking model.
# It uses the parameters and models of the environment.csv file within this folder

# TODO: Maybe we do need to check the FE in the gene region.
## :
# 1.- Load csv files as dict
# 2.- Define params to calibrate and ranges (for the moment, the distance and the foldchange)
# 3.- Similar to calibrate_topos.py, add the new params to the dict to calibrate the enzyme.
# 4.- Define systems to vary. For the moment, don't vary anything. You need to find a way to define the circuit
#     by hand. That means that you don't need a csv. Maybe you need to create a new function
#     So for the moment, don't vary anything in the system, but keep in mind that you should vary
#     Domain size (small, medium, large) and promoter strength (weak, medium, strong)
#    O al chile, you can just three different input files at the moment, instead of wasting time coding additional
#    stuff. O tal vez no es tan dificil definir un circuit, tal vez solo darle un dict (como csv).
# 5.- Write function that performs the calibration process. Needs to run simulations. Calculate counts.
#     Fit distribution. Calculate avg Foldenrichment, and cross correlation.
#     ** Note that you need a curve to compare the fold enrichment. So for each domain size, you'll have a different
#        density.
# 6.- You can later make everything a list (reference and everything) so you can vary the system variables

# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 60 #12#24  #24#8 #96 # 120
tests = 60  # number of tests for parametrization

# Simulation conditions
# --------------------------------------------------------------
file_out = 'calibration_RNAPTracking'
dt = 0.25
initial_time = 0
final_time = 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Reference - It is the reference density of topos when there is no gene that we will use to calculate the
#             fold enrichment.
# --------------------------------------------------------------
reference_path = '../noRNAP/'  # Path with the reference kdes

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = '../sites.csv'
enzymes_filename = None #'../enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'topoIRNAPtrack-uniform'
series = True
continuation = False

initial_sigma = -0.044  # The initial superhelical density
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
target_FE = 1.68  # Target fold-enrichment
target_CO = 1.0  # Target correlation between topo I and RNAP densities.
x_spacing = 10.0 # The spacing I want at the moment of doing interpolation.

# nbins is the number of bins to use when calculating the kde
target_dict = {'target_FE': target_FE, 'target_CO': target_CO, 'target_gene': 'reporter',
               'enzymes_names': ['RNAP', 'topoI', 'gyrase']}

# ----------------------------------------------------------------------------------------------------------------------
# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------

def objective_function(params):
    # We need to prepare the inputs.
    # At the moment, we only have one system.

    # Global dictionaries
    # ------------------------------------------
    global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                   'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                   'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                   'frames': frames, 'dt': dt, 'n_simulations': n_simulations, # 'initial_sigma': initial_sigma,
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
    my_objective, output_dict = tct.single_case_RNAPTracking_calibration(global_dict_list, variations_list,
                                                                         list_reference, n_simulations,
                                                                         target_dict)

    return my_objective


# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------

# Let's load the circuit, so we can extract some information
# -----------------------------------
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

# Get target site
target_gene = [site for site in my_circuit.site_list if site.name == target_dict['target_gene']][0]

# Define x-axes
x_system = tct.get_interpolated_x(1, my_circuit.size)
x_gene = tct.get_interpolated_x(target_gene.start, target_gene.end)

reference_dict = {}
for name in target_dict['enzymes_names']:

    # Load reference file
    if name != 'RNAP':
        kde_ref = np.loadtxt(reference_path + 'reference_'+ name + '.txt')

        # And get interpolated data - I do it here, so we don't have to do it again in the parallelization and save some
        # time. It is saved to the reference_dict by the way
        reference_dict[name] = tct.get_interpolated_kde(kde_ref[:, 0], kde_ref[:, 1], x_system)

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
    sys.stdout = f##

    # Your code that prints to the screen
    print("Hello, this is the info file for the calibration of Topo I Tracking RNAP.")
    print("Topo I Binding Model = " + topoI_binding_model_name)
    print('Ran ' + str(n_simulations) + ' simulations per test. ')
    print('Number of tests = ' + str(tests))

    best = fmin(
        fn=objective_function,  # Objective Function to optimize
        space=space,  # Hyperparameter's Search Space
        algo=tpe.suggest,  # Optimization algorithm (representative TPE)
        max_evals=tests,  # Number of optimization attempts
        trials = trials

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

    tdi = trials.trials[n] # dictionary with results for test n
    lo = trials.trials[n]['result']['loss'] # loss
    va = trials.trials[n]['misc']['vals'] #values
    # Add a new row using append method
    new_row = pd.DataFrame({
        'test': n, 'loss': lo, 'RNAP_dist': va['RNAP_dist'], 'fold_change': va['fold_change'],
        'k_on': k_on, 'width': width, 'threshold': threshold
    })
#    params_df.append(new_row, ignore_index=True)
    params_df = pd.concat([params_df, new_row], ignore_index=True)


params_df.to_csv('values.csv', index=False, sep=',')
