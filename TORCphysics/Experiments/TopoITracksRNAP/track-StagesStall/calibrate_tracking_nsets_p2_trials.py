import sys
sys.path.append("/users/ph1vvb")
import numpy as np
import pandas as pd
from hyperopt import tpe, hp, fmin, Trials
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit, params
import pickle

# Description
# --------------------------------------------------------------
# Following same logic than calibrate_tracking_nsets, this script runs the calibration using the function:
# single_case_RNAPTracking_calibration_nsets_p2()
# This function uses a different parallelization strategy to speed up the calibration process.
# It does this by doing two parallelization pools.

# This script launches the calibration process for the TopoIRNAPTracking binding model, SpacerBinding binding model,
# RNAPStagesStall effect model and the RNAPStagesSimpleUnbinding unbinding model.
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
#n_workers = 16 # Total number of workers (cpus) - For us testing in this machine
n_workers = 64 # For stanage
#n_sets = 1 #7  # Number of outer sets
n_sets = 7 # - For stange
n_subsets = 4   # Number of simulations per set - One KDE for each n_subset within the n_set
                  # As it says upthere, there are in total n_subsets * n_sets number of KDEs per test
                  # A number of n_inner_workers parallel simulations are launch to calculate each individual KDE.
n_inner_workers = n_workers // (n_sets+1)  # Number of workers per inner pool
                                           # +1 because one worker is spent in creating the outer pool
                                           # I think this number of n_inner_workers is the number of parallel simulations
                                           # ran that are used to calculate each individual KDE.
#tests = 4 #400  # number of tests for parametrization
tests = NTESTS # number of tests for parametrization
dt = NDT
final_time = NTIME

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
initial_time = 0
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'SIMNUMBER-calibration_avg-RNAPTracking_nsets_p2_small_dt'+str(dt)+'_smalldist'

# Reference - It is the reference density of topos when there is no gene that we will use to calculate the
#             fold enrichment.
# --------------------------------------------------------------
#reference_path = '../noRNAP/'  # Path with the reference kdes
reference_path = '../noRNAP/avg_'  # Path with the reference kdes
reference_RNAP = '../RNAP_TU.txt'

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = 'sites.csv' # '../sites.csv'
enzymes_filename = None  #'../enzymes.csv'
#environment_filename = 'environment.csv'
#environment_filename = 'environment_dt'+str(dt)+'.csv'
environment_filename = 'environment_avg_dt'+str(dt)+'.csv'
output_prefix = 'topoIRNAPtrack-StagesStall_dt'+str(dt)
series = True
continuation = False

# Let's asume it starts at a value where topos are equilibrated, so we assume the steady state.

# Models to calibrate
# -----------------------------------
# Topoisomerase I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'TopoIRecognitionRNAPTracking'
# These are parameters for the Recognition curve
#topoI_params_csv = '../calibration_dt'+str(dt)+'_topoI.csv'
topoI_params_csv = '../topoI_rec_avg_dt'+str(dt)+'.csv'

# Site
gene_name = 'reporter'
gene_type = 'site'
gene_binding_model_name = 'SpacerBinding'
spacer_oparams = {'superhelical_op': -0.06, 'spacer': 17}

# RNAP
RNAP_name = 'RNAP'
RNAP_type = 'environmental'
RNAP_effect_model_name = 'RNAPStagesStall'
RNAP_unbinding_model_name = 'RNAPStagesSimpleUnbinding'
RNAP_oparams = {'width': 0.005, 'threshold': -0.042, # Sam Meyer's PROMOTER CURVE (parameters taken from Houdaigi NAR 2019)
                'velocity': params.v0, 'kappa': params.RNAP_kappa, 'stall_torque': params.stall_torque}

# RANGES FOR RANDOM SEARCH
# -----------------------------------
# Topo model
RNAP_dist_min = 20
RNAP_dist_max = 120 #500
fold_change_min = .1
fold_change_max = 30#50

# Reporter
spacer_kon_min = 0.03
spacer_kon_max = .3

# RNAP - rates
k_closed_min = 0.03
k_closed_max = 0.3
k_open_min = 0.03
k_open_max = 0.3
k_ini_min = 0.05
k_ini_max = 0.3
gamma_min = 0.03
gamma_max = 1.0
k_off_min = 0.03
k_off_max = 0.3


# TARGETS FOR OPTIMIZATION
# -----------------------------------
target_FE = 1.238 # 3rd System#1.38 #1.68  # Target fold-enrichment
target_CO = 0.944 # 3rd System #1.0  # Target correlation between topo I and RNAP densities.
target_RNAP_CO = 1.0 # target correlation I want between RNAP KDE from simulations and ChIP-Seq
x_spacing = 10.0  # The spacing I want at the moment of doing interpolation.

# nbins is the number of bins to use when calculating the kde
target_dict = {'target_FE': target_FE, 'target_CO': target_CO, 'target_RNAP_CO': target_RNAP_CO,
               'target_gene': 'reporter',
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

    # Gene (site)
    name = gene_name
    object_type = gene_type
    binding_model_name = gene_binding_model_name
    binding_oparams = {'k_on': params['spacer_kon'], 'spacer': spacer_oparams['spacer'],
                       'superhelical_op': spacer_oparams['superhelical_op'],}

    gene_variation = {'name': name, 'object_type': object_type,
                      'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams}

    # RNAP (environmental - effect model + unbinding model)
    name = RNAP_name
    object_type = RNAP_type
    effect_model_name = RNAP_effect_model_name
    effect_oparams = { 'k_closed':params['k_closed'], 'k_open':params['k_open'],
                       'width': RNAP_oparams['width'], 'threshold': RNAP_oparams['threshold'],
                       'k_ini': params['k_ini'],
                       'velocity': RNAP_oparams['velocity'], 'gamma': params['gamma'],
                       'stall_torque': RNAP_oparams['stall_torque'], 'kappa': RNAP_oparams['kappa']}
    unbinding_model_name = RNAP_unbinding_model_name
    unbinding_oparams = {'k_off': params['k_off']}
    RNAP_variation = {'name': name, 'object_type': object_type,
                      'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                      'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams}

    # Create lists of conditions for each system
    # ------------------------------------------

    # Global dictionaries
    global_dict_list = [global_dict]

    # List of lists of variations - It is a list of lists because we have a list for each experiment.
    #                               So for 1 experiment and three variations, the list is on the form [[v1,v2,v3]],
    #                               if it were 2 experiments and 3 variations each: [[v1.1,v1.2,v1.3], [v2.1,v2.2,v2.3]]
    variations_list = [[topoI_variation, gene_variation, RNAP_variation]]

    # Arrays with position densities to calculate fold change
    list_reference = [reference_dict]

    # Finally, run objective function.
    # ------------------------------------------
    my_objective, output_dict = tct.single_case_RNAPTracking_calibration_nsets_2scheme_plus_RNAP_odict(global_dict_list,
                                                                                                 variations_list,
                                                                                                 list_reference,
                                                                                                 target_dict)
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
        kde_ref = np.loadtxt(reference_path + 'reference_' + name + '_dt'+str(dt)+'.txt')

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
    'fold_change': hp.uniform('fold_change', fold_change_min, fold_change_max),

    # RNAP params
    'spacer_kon': hp.uniform('spacer_kon', spacer_kon_min, spacer_kon_max),

    'k_closed': hp.uniform('k_closed', k_closed_min, k_closed_max),
    'k_open': hp.uniform('k_open', k_open_min, k_open_max),
    'k_ini': hp.uniform('k_ini', k_ini_min, k_ini_max),
    'gamma': hp.uniform('gamma', gamma_min, gamma_max),

    'k_off': hp.uniform('k_off', k_off_min, k_off_max)
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

# Let's save the trials object
# --------------------------------------------------------------------------
with open(file_out+'-trials.pkl', 'wb') as file:
    pickle.dump(trials, file)

# Save model!
# --------------------------------------------------------------------------
# Topo I
k_on = float(topoI_params['k_on'][0])
width = float(topoI_params['width'][0])
threshold = float(topoI_params['threshold'][0])

# RNAP-spacer
#superhelical_op = spacer_oparams['superhelical_op']
#spacer = spacer_oparams['spacer']
#RNAP_width = RNAP_oparams['width']
#RNAP_threshold = RNAP_oparams['threshold']
#kappa = RNAP_oparams['kappa']
#velocity = RNAP_oparams['velocity']
#stall_torque = RNAP_oparams['stall_torque']

# This one will have all the values of the RNAPTracking Model for the enzyme calibration
# Topo I
name = 'topoI-'
complete_df = pd.DataFrame(columns=['k_on', 'width', 'threshold', 'RNAP_dist', 'fold_change'])
complete_df['RNAP_dist'] = best_df['RNAP_dist']
complete_df['fold_change'] = best_df['fold_change']
complete_df['k_on'] = k_on
complete_df['width'] = width
complete_df['threshold'] = threshold
complete_df['k_off'] = float(topoI_params['k_off'][0])
complete_df['k_cat'] = float(topoI_params['k_cat'][0])

complete_df.to_csv(name+file_out + '.csv', index=False, sep=',')

# RNAP
name = 'RNAP-'
complete_df = pd.DataFrame(columns=['k_on', 'superhelical_op', 'spacer',
                                    'k_closed', 'k_open', 'width', 'threshold',
                                    'k_ini', 'gamma', 'kappa', 'velocity', 'stall_torque',
                                    'k_off'])
complete_df['k_on'] = best_df['spacer_kon']
complete_df['superhelical_op'] = spacer_oparams['superhelical_op']
complete_df['spacer'] = spacer_oparams['spacer']
complete_df['k_closed'] = best_df['k_closed']
complete_df['k_open'] = best_df['k_open']
complete_df['width'] = RNAP_oparams['width']
complete_df['threshold'] = RNAP_oparams['threshold']
complete_df['k_ini'] = best_df['k_ini']
complete_df['gamma'] = best_df['gamma']
complete_df['kappa'] = RNAP_oparams['kappa']
complete_df['velocity'] = RNAP_oparams['velocity']
complete_df['stall_torque'] = RNAP_oparams['stall_torque']
complete_df['k_off'] = best_df['k_off']
complete_df.to_csv(name+file_out + '.csv', index=False, sep=',')

# Let's save trials info (params and loses)
# --------------------------------------------------------------------------
# Topo I
name = 'topoI-'
params_df = pd.DataFrame(columns=['test', 'loss', 'RNAP_dist', 'fold_change'])
for n in range(tests):
    tdi = trials.trials[n]  # dictionary with results for test n
    lo = trials.trials[n]['result']['loss']  # loss
    va = trials.trials[n]['misc']['vals']  #values
    # Add a new row using append method
    new_row = pd.DataFrame({
        'test': n, 'loss': lo, 'RNAP_dist': va['RNAP_dist'], 'fold_change': va['fold_change'],
        'k_on': k_on, 'width': width, 'threshold': threshold, 'k_off': topoI_params['k_off'][0],
        'k_cat': topoI_params['k_cat'][0]
    })
    #    params_df.append(new_row, ignore_index=True)
    params_df = pd.concat([params_df, new_row], ignore_index=True)

params_df.to_csv(name+file_out+'-values.csv', index=False, sep=',')

# RNAP
name = 'RNAP-'
params_df = pd.DataFrame(columns=['test', 'loss'])
for n in range(tests):
    tdi = trials.trials[n]  # dictionary with results for test n
    lo = trials.trials[n]['result']['loss']  # loss
    va = trials.trials[n]['misc']['vals']  #values
    # Add a new row using append method
    new_row = pd.DataFrame({
        'test': n, 'loss': lo,
        'k_on': va['spacer_kon'], 'superhelical_op': spacer_oparams['superhelical_op'], 'spacer': spacer_oparams['spacer'],
        'k_closed': va['k_closed'], 'k_open': va['k_open'], 'width': RNAP_oparams['width'], 'threshold': RNAP_oparams['threshold'],
        'k_ini': va['k_ini'], 'gamma': va['gamma'], 'kappa': RNAP_oparams['kappa'], 'velocity': RNAP_oparams['velocity'],
        'stall_torque': RNAP_oparams['stall_torque'], 'k_off': va['k_off']
    })

    params_df = pd.concat([params_df, new_row], ignore_index=True)

params_df.to_csv(name+file_out+'-values.csv', index=False, sep=',')

# Let's run the function once more with the best params to produce the data so we can then just plot it.
# --------------------------------------------------------------------------
objective, output_dict = objective_function(params=best, calibrating=False)

# Save the dictionary to a file
with open(file_out+'.pkl', 'wb') as file:
    pickle.dump(output_dict, file)
