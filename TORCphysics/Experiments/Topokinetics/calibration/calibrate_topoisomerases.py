#import sys
#sys.path.append("/users/ph1vvb")  # For stanage
import numpy as np
from hyperopt import tpe, hp, fmin, Trials
import pandas as pd
from TORCphysics import topo_calibration_tools as tct

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This is just a test to reproduce the global supercoiling response curves from the paper:
# Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching
# This one controlls the fast unbinding - checks that the k_off rate is smaller than the k_on in general,
# Because we have to control the number of molecules through the simulation, if not, eventually everything will
# be full of topoisomerases.

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0
tests = 10000   # number of tests for parametrization

# For parallelization and calibration
n_simulations = 100 # For stanage

initial_time = 0
# Let's do it for 400s to add more weight to the curve and not the plateau
final_time = 500  #500 # 600
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'rep-calibration_dt' + str(dt)

# For the simulation
circuit_filename = '../circuit.csv'
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
environment_filename = '../environment.csv'

# Concentrations in nM
DNA_concentration = 0.75
gyrase_concentration = 44.6
topoI_concentration = 17.0

# MM kinetics
K_M_topoI = 1.5
k_cat_topoI = .0023
v_max_topoI = k_cat_topoI * topoI_concentration
K_M_gyrase = 2.7
k_cat_gyrase = .0011
v_max_gyrase = k_cat_gyrase * gyrase_concentration

# Superhelical values (sigma) for each case
sigma_0_topo = -0.11  #-0.076  # Approximately -20 supercoils according the paper
sigma_0_gyrase = 0.0  # We suppose this one.
sigma_f_gyrase = -0.11  # We also assume this one, which is the maximum at which gyrase acts.
# At this value the torque is too strong.

output_prefix = 'test0'
series = True
continuation = False


# Models to calibrate to calibrate
# -----------------------------------
# Topoisomerase I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'TopoIRecognition'
topoI_effect_model_name = 'TopoILinear'
topoI_unbinding_model_name = 'PoissonUnBinding'

# Gyrase
gyrase_name = 'gyrase'
gyrase_type = 'environmental'
gyrase_binding_model_name = 'GyraseRecognition'
gyrase_effect_model_name = 'GyraseLinear'
gyrase_unbinding_model_name = 'PoissonUnBinding'

# RANGES FOR RANDOM SEARCH
# -----------------------------------
# TopoI ranges
# I think this means that, in the max k_on rate, 1 topo binds every second, and in the minimum, one topo binds every
# ~ 100 seconds; this for the whole circuit. - This for the small environment (realistic).
# The true k_on of the whole circuit depends on the size of the plasmid, the size of the enzyme, and the concentration.
#k_on_min_topoI = 0.0000042 #0.001
k_on_min_topoI = 0.000042 #0.001
#k_on_max_topoI = 4*0.00042 #0.1
k_on_max_topoI = 0.00042 #0.1
k_off_min_topoI = 0.1
k_off_max_topoI = 1.0
k_cat_min_topoI = 0.0  # Ranges to vary k_cat
k_cat_max_topoI = 20.0
# topoI new ranges
#thres max=0.03, min=-0.07
#width max=0.0125, min=0.00125
width_min_topoI = 0.00125 #0.001
width_max_topoI = 0.0125 #0.3  #0.05
threshold_min_topoI = -0.07#  -0.3  # 0.001  #-0.05
threshold_max_topoI = 0.03 #0.3  #0.05  #-0.001

# Gyrase ranges
# For gyrase, the k_on of gyrase for the whole circuit is similar to the topoI, which goes from 1 bind every second
# to 1 event every 100 seconds
#k_on_min_gyrase = 0.0000024#0.001
k_on_min_gyrase = 0.000024#0.001
#k_on_max_gyrase = 4*0.00024#0.1
k_on_max_gyrase = 0.00024#0.1
k_off_min_gyrase = 0.1
k_off_max_gyrase = 1.0
k_cat_min_gyrase = 0.0  # Ranges to vary k_cat
k_cat_max_gyrase = 20.0
sigma0_min_gyrase = -.15  #-0.3
sigma0_max_gyrase = -0.07  #0.0
# gyrase new ranges
#thres max=0.04, min=-0.06
#width max=0.0125, min=0.00125
width_min_gyrase = 0.00125 #0.001
width_max_gyrase = 0.0125 #0.3  #.05  # 1.0
threshold_min_gyrase = -0.06 #-0.3  #0.0 #0.001
threshold_max_gyrase = 0.04 #0.3  # 0.1  # 1.0

# Sizes
plasmid_size = 2757
topoI_size = 20
gyrase_size = 30
topoI_ngrid = int(plasmid_size/topoI_size)
gyrase_ngrid = int(plasmid_size/gyrase_size)

# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the mean superhelical
# density for each substrate concentration

def objective_function(params):

    # First, let's test the k_on/k_off ratio
    # ------------------------------------------------------------------
    test_topoI_kon = (params['k_on_topoI'] * topoI_ngrid * topoI_concentration)
    test_gyrase_kon =  (params['k_on_gyrase'] * gyrase_ngrid * gyrase_concentration)

    topoI_ratio = test_topoI_kon / params['k_off_topoI']
    gyrase_ratio = test_gyrase_kon / params['k_off_gyrase']

    # Test ratio - if it's bigger than one, then we don't want it
    if topoI_ratio > 1.0 or gyrase_ratio > 1.0:
        return 100 # Big objective

    # We need to prepare the inputs.
    # This time we have three different systems:
    # 1.- Topoisomerase I acting on supercoiled DNA
    # 2.- Gyrase acting on Relaxed DNA
    # 3.- Both enzymes acting on supercoiled/Relaxed DNA

    # Global dictionaries
    # ------------------------------------------
    global_dict_topoI = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                         'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                         'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                         'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma_0_topo,
                         'DNA_concentration': 0.0}

    global_dict_gyrase = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                          'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                          'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                          'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma_0_gyrase,
                          'DNA_concentration': 0.0}

    global_dict_both_sc = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                           'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                           'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                           'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma_0_topo,
                           'DNA_concentration': 0.0}

    global_dict_both_rx = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                           'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                           'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                           'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma_0_gyrase,
                           'DNA_concentration': 0.0}

    # Variation dictionaries
    # ------------------------------------------

    # Topoisomerase I
    name = topoI_name
    object_type = topoI_type
    binding_model_name = topoI_binding_model_name
    binding_oparams = {'k_on': params['k_on_topoI'], 'width': params['width_topoI'],
                       'threshold': params['threshold_topoI']}
    effect_model_name = topoI_effect_model_name
    effect_oparams = {'k_cat': params['k_cat_topoI']}
    unbinding_model_name = topoI_unbinding_model_name
    unbinding_oparams = {'k_off': params['k_off_topoI']}
    concentration = topoI_concentration  # / mol_concentration  # Because this is the reference.

    topoI_variation = {'name': name, 'object_type': object_type,
                       'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                       'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                       'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                       'concentration': concentration}

    # Gyrase
    name = gyrase_name
    object_type = gyrase_type
    binding_model_name = gyrase_binding_model_name
    binding_oparams = {'k_on': params['k_on_gyrase'], 'width': params['width_gyrase'],
                       'threshold': params['threshold_gyrase']}
    effect_model_name = gyrase_effect_model_name
    effect_oparams = {'k_cat': params['k_cat_gyrase'], 'sigma0': params['sigma0_gyrase']}
    unbinding_model_name = gyrase_unbinding_model_name
    unbinding_oparams = {'k_off': params['k_off_gyrase']}
    concentration = gyrase_concentration  # / mol_concentration  # Because this is the reference.

    gyrase_variation = {'name': name, 'object_type': object_type,
                        'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                        'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                        'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                        'concentration': concentration}

    # Create lists of conditions for each system
    # ------------------------------------------

    # Global dictionaries
    global_dict_list = [global_dict_topoI, global_dict_gyrase, global_dict_both_sc, global_dict_both_rx]

    # List of lists of variations
    variations_list = [[topoI_variation], [gyrase_variation],
                       [topoI_variation, gyrase_variation],
                       [topoI_variation, gyrase_variation]
                       ]

    # Arrays with global superhelical densities
    list_sigmas = [topoI_sigma, gyrase_sigma, both_sigma_sc, both_sigma_rx]

    # Finally, run objective function. run_objective_function will process our conditions
    # ------------------------------------------
    my_objective, simulation_superhelicals = tct.run_objective_function(global_dict_list=global_dict_list,
                                                                        variations_list=variations_list,
                                                                        exp_superhelicals=list_sigmas,
                                                                        n_simulations=n_simulations)
    return my_objective


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------
# Build experimental curve for TOPO I
# -----------------------------------------------------
# Kinetics: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI
# Product = Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs; which initially is the same as the DNA concentration

# Integrate MM kinetics
# ------------------------------------------
# Initially, there's no relaxed DNA, and all the supercoiled DNA concentration corresponds to the plasmid conc.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_topoI(vmax=v_max_topoI, KM=K_M_topoI,
                                                      Supercoiled_0=DNA_concentration, Relaxed_0=0.0,
                                                      frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
# Note that this sigmaf is not at the end of the simulation, but the sigma at which there is 0 Relaxed DNA.
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
topoI_sigma = sigma

# -----------------------------------------------------
# Build experimental curve for Gyrase
# -----------------------------------------------------
# Kinetics: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Product = Supercoiled DNA
# Substrate = Relaxed DNA; which initially is the same as the DNA concentration

# Integrate MM kinetics
# ------------------------------------------
# Initially, there's no supercoiled DNA, and all of the relaxed DNA concentration corresponds
# to the plasmid concentration.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_gyrase(vmax=v_max_gyrase, KM=K_M_gyrase,
                                                       Supercoiled_0=0.0, Relaxed_0=DNA_concentration,
                                                       frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
gyrase_sigma = sigma

# -----------------------------------------------------
# Build experimental curve for system with both Topo I and Gyrase (from Sc state)
# -----------------------------------------------------
# Kinetics Gyrase: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Kinetics Topoisomerase: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI

# Integrate MM kinetics
# ------------------------------------------
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_both_T_G(vmax_topoI=v_max_topoI, vmax_gyrase=v_max_gyrase,
                                                         KM_topoI=K_M_topoI, KM_gyrase=K_M_gyrase,
                                                         Supercoiled_0=DNA_concentration, Relaxed_0=0.0,
                                                         frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
both_sigma_sc = sigma

# -----------------------------------------------------
# Build experimental curve for system with both Topo I and Gyrase (from Rx state)
# -----------------------------------------------------
# Integrate MM kinetics
# ------------------------------------------
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_both_T_G(vmax_topoI=v_max_topoI, vmax_gyrase=v_max_gyrase,
                                                         KM_topoI=K_M_topoI, KM_gyrase=K_M_gyrase,
                                                         Supercoiled_0=0.0, Relaxed_0=DNA_concentration,
                                                         frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
both_sigma_rx = sigma

# Optimization
# -----------------------------------------------------
trials = Trials()

space = {

    # Topo I params
    'k_cat_topoI': hp.uniform('k_cat_topoI', k_cat_min_topoI, k_cat_max_topoI),
    'k_on_topoI': hp.uniform('k_on_topoI', k_on_min_topoI, k_on_max_topoI),
    'k_off_topoI': hp.uniform('k_off_topoI', k_off_min_topoI, k_off_max_topoI),
    'width_topoI': hp.uniform('width_topoI', width_min_topoI, width_max_topoI),
    'threshold_topoI': hp.uniform('threshold_topoI', threshold_min_topoI, threshold_max_topoI),

    # Gyrase params
    'k_cat_gyrase': hp.uniform('k_cat_gyrase', k_cat_min_gyrase, k_cat_max_gyrase),
    'k_on_gyrase': hp.uniform('k_on_gyrase', k_on_min_gyrase, k_on_max_gyrase),
    'k_off_gyrase': hp.uniform('k_off_gyrase', k_off_min_gyrase, k_off_max_gyrase),
    'sigma0_gyrase': hp.uniform('sigma0_gyrase', sigma0_min_gyrase, sigma0_max_gyrase),
    'width_gyrase': hp.uniform('width_gyrase', width_min_gyrase, width_max_gyrase),
    'threshold_gyrase': hp.uniform('threshold_gyrase', threshold_min_gyrase, threshold_max_gyrase)
}

# Save the current standard output
#original_stdout = sys.stdout

# Define the file where you want to save the output
output_file_path = file_out + '.info'

# Open the file in write mode
with open(output_file_path, 'w') as f:
    # Redirect the standard output to the file
#    sys.stdout = f

    # Your code that prints to the screen
    print("Hello, this is the info file for the calibration of Topo I and Gyrase Models.")
    print("Launching calibration for dt="+str(dt))
    print("Topo I Binding Model = " + topoI_binding_model_name)
    print("Topo I Effect Model = " + topoI_effect_model_name)
    print("Topo I Unbinding Model = " + topoI_unbinding_model_name)
    print("Gyrase Binding Model = " + gyrase_binding_model_name)
    print("Gyrase Effect Model = " + gyrase_effect_model_name)
    print("Gyrase Unbinding Model = " + gyrase_unbinding_model_name)
    print('Ran ' + str(n_simulations) + ' simulations per test. ')
    print('Number of tests = ' + str(tests))

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
best_df.to_csv(file_out + '.csv', index=False, sep=',')

# Let's save trials info (params and loses)
# --------------------------------------------------------------------------
params_df = pd.DataFrame(columns=['test', 'loss',
                                  'k_on_topoI', 'k_off_topoI', 'k_cat_topoI',
                                  'width_topoI', 'threshold_topoI',

                                  'k_on_gyrase', 'k_off_gyrase', 'k_cat_gyrase',
                                  'width_gyrase', 'threshold_gyrase', 'sigma0_gyrase'])

for n in range(tests):
    tdi = trials.trials[n]  # dictionary with results for test n
    lo = trials.trials[n]['result']['loss']  # loss
    va = trials.trials[n]['misc']['vals']  #values
    # Add a new row using append method
    new_row = pd.DataFrame({
        'test': n, 'loss': lo,
        'k_on_topoI': va['k_on_topoI'], 'k_off_topoI': va['k_off_topoI'], 'k_cat_topoI': va['k_cat_topoI'],
        'width_topoI': va['width_topoI'], 'threshold_topoI': va['threshold_topoI'],
        'k_on_gyrase': va['k_on_gyrase'], 'k_off_gyrase': va['k_off_gyrase'], 'k_cat_gyrase': va['k_cat_gyrase'],
        'width_gyrase': va['width_gyrase'], 'threshold_gyrase': va['threshold_gyrase'],
        'sigma0_gyrase': va['sigma0_gyrase']
    })
    #    params_df.append(new_row, ignore_index=True)
    params_df = pd.concat([params_df, new_row], ignore_index=True)

params_df.to_csv(file_out+'_values.csv', index=False, sep=',')

# Let's save it for each enzyme
topo_df = pd.DataFrame(columns=['k_on', 'k_off', 'k_cat', 'width', 'threshold'])
topo_df['k_on'] = best_df['k_on_topoI']
topo_df['k_off'] = best_df['k_off_topoI']
topo_df['k_cat'] = best_df['k_cat_topoI']
topo_df['width'] = best_df['width_topoI']
topo_df['threshold'] = best_df['threshold_topoI']
topo_df.to_csv(file_out+'_topoI.csv', index=False, sep=',')

gyrase_df = pd.DataFrame(columns=['k_on', 'k_off', 'k_cat', 'width', 'threshold', 'sigma0'])
gyrase_df['k_on'] = best_df['k_on_gyrase']
gyrase_df['k_off'] = best_df['k_off_gyrase']
gyrase_df['k_cat'] = best_df['k_cat_gyrase']
gyrase_df['width'] = best_df['width_gyrase']
gyrase_df['threshold'] = best_df['threshold_gyrase']
gyrase_df['sigma0'] = best_df['sigma0_gyrase']
gyrase_df.to_csv(file_out+'_gyrase.csv', index=False, sep=',')
