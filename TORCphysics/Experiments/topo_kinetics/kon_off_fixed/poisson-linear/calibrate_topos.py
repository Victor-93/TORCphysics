import sys
import numpy as np
from hyperopt import tpe, hp, fmin
import pandas as pd
from TORCphysics import topo_calibration_tools as tct

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This is to reproduce the global supercoiling response curves from the paper:
# Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching

# This time will deduce k_on & k_off based on the kinetic measurements and assuming that k_off = 0.5 (2 seconds),
# based on:
# Single-molecule imaging of DNA gyrase activity in living Escherichia coli

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 0.25
initial_time = 0
# Let's do it for 400s to add more weight to the curve and not the plateau
final_time = 400  #500 # 600
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'calibration'

# For the simulation
circuit_filename = '../../circuit.csv'
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
environment_filename = '../../environment.csv'

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

# For parallelization and calibration
n_simulations = 96  #60#84  # 60 #48 #120
tests = 2000#10  #400#100 #400  # 10  # 100  # number of tests for parametrization

# Models to calibrate to calibrate
# -----------------------------------
# Topoisomerase I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'PoissonBinding'
topoI_effect_model_name = 'TopoILinear'
topoI_unbinding_model_name = 'PoissonUnBinding'

# Gyrase
gyrase_name = 'gyrase'
gyrase_type = 'environmental'
gyrase_binding_model_name = 'PoissonBinding'
gyrase_effect_model_name = 'GyraseLinear'
gyrase_unbinding_model_name = 'PoissonUnBinding'

# SOME SIZES
# -----------------------------------
plasmid_length = 2757 # Make sure is the same as in circuit.csv
topoI_size = 60
gyrase_size = 160
n_grid_topoI = int(plasmid_length/topoI_size)
n_grid_gyrase = int(plasmid_length/gyrase_size)

# FIXED VALUES
# -----------------------------------
k_off_topoI = 0.5
k_off_gyrase = 0.5

# Deduce k_on based on k_off
k_on_topoI = ((k_off_topoI + k_cat_topoI) / K_M_topoI)/(n_grid_topoI*topoI_concentration)
k_on_gyrase = ((k_off_gyrase + k_cat_gyrase) / K_M_gyrase)/(n_grid_gyrase*gyrase_concentration)

# RANGES FOR RANDOM SEARCH
# -----------------------------------
# TopoI ranges
k_cat_min_topoI = 1.0  # Ranges to vary k_cat
k_cat_max_topoI = 20.0

# Gyrase ranges
k_cat_min_gyrase = 1.0  # Ranges to vary k_cat
k_cat_max_gyrase = 20.0

# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the mean superhelical
# density for each substrate concentration

def objective_function(params):
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

    global_dict_both = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                        'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                        'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                        'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma_0_topo,
                        'DNA_concentration': 0.0}

    # Variation dictionaries
    # ------------------------------------------

    # Topoisomerase I
    name = topoI_name
    object_type = topoI_type
    binding_model_name = topoI_binding_model_name
    binding_oparams = {'k_on':k_on_topoI}
    effect_model_name = topoI_effect_model_name
    effect_oparams = {'k_cat': params['k_cat_topoI']}
    unbinding_model_name = topoI_unbinding_model_name
    unbinding_oparams = {'k_off': k_off_topoI}
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
    binding_oparams = {'k_on': k_on_gyrase}
    effect_model_name = gyrase_effect_model_name
    effect_oparams = {'k_cat': params['k_cat_gyrase'], 'sigma0': sigma_f_gyrase}  #params['sigma0_gyrase']}
    unbinding_model_name = gyrase_unbinding_model_name
    unbinding_oparams = {'k_off': k_off_gyrase}
    concentration = gyrase_concentration  # / mol_concentration  # Because this is the reference.

    gyrase_variation = {'name': name, 'object_type': object_type,
                        'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                        'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                        'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                        'concentration': concentration}

    # Create lists of conditions for each system
    # ------------------------------------------

    # Global dictionaries
    global_dict_list = [global_dict_topoI, global_dict_gyrase, global_dict_both]

    # List of lists of variations
    variations_list = [[topoI_variation], [gyrase_variation], [topoI_variation, gyrase_variation]]

    # Arrays with global superhelical densities
    list_sigmas = [topoI_sigma, gyrase_sigma, both_sigma]

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
#gyrase_sigma = tct.gyrase_to_sigma(Relaxed=relaxed_DNA, DNA_concentration=DNA_concentration,
#                                   sigma0=sigma_0_gyrase, sigmaf=sigma_f_gyrase)

# -----------------------------------------------------
# Build experimental curve for system with both Topo I and Gyrase
# -----------------------------------------------------
# Kinetics Gyrase: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Kinetics Topoisomerase: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI

# Integrate MM kinetics
# ------------------------------------------

# Initially, there's no supercoiled DNA, and all of the relaxed DNA concentration corresponds
# to the plasmid concentration.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_both_T_G(vmax_topoI=v_max_topoI, vmax_gyrase=v_max_gyrase,
                                                         KM_topoI=K_M_topoI, KM_gyrase=K_M_gyrase,
                                                         Supercoiled_0=DNA_concentration, Relaxed_0=0.0,
                                                         frames=frames, dt=dt)
ratio = relaxed_DNA[-1] / DNA_concentration
sigmaf = sigma_0_topo * ratio
# Translate to superhelical density
# ------------------------------------------
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
both_sigma = sigma

#both_sigma = tct.both_T_G_to_sigma(Relaxed=relaxed_DNA, Relaxed_final=relaxed_DNA[-1],
#                                   sigma0=sigma_0_topo, sigmaf=sigmaf)

# Optimization
# -----------------------------------------------------
space = {

    # Topo I params
    'k_cat_topoI': hp.uniform('k_cat_topoI', k_cat_min_topoI, k_cat_max_topoI),

    # Gyrase params
    'k_cat_gyrase': hp.uniform('k_cat_gyrase', k_cat_min_gyrase, k_cat_max_gyrase)
}

# Save the current standard output
original_stdout = sys.stdout
# Define the file where you want to save the output
output_file_path = file_out + '.info'

# Open the file in write mode
with open(output_file_path, 'w') as f:
    # Redirect the standard output to the file
    sys.stdout = f

    # Your code that prints to the screen
    print("Hello, this is the info file for the calibration of Topo I and Gyrase Models.")
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
        max_evals=tests  # Number of optimization attempts
    )

    print(" ")
    print("Optimal parameters found from random search: ")
    print(best)

best_df = pd.DataFrame.from_dict([best])

# Let's save it for each enzyme
topo_df = pd.DataFrame(columns=['k_on', 'k_off', 'k_cat'])
topo_df['k_cat'] = best_df['k_cat_topoI']
topo_df['k_on'] = k_on_topoI
topo_df['k_off'] = k_off_topoI
topo_df.to_csv('calibration_topoI.csv', index=False, sep=',')

gyrase_df = pd.DataFrame(columns=['k_on', 'k_off', 'k_cat', 'sigma0'])
gyrase_df['k_cat'] = best_df['k_cat_gyrase']
gyrase_df['sigma0'] = sigma_f_gyrase
gyrase_df['k_on'] = k_on_gyrase
gyrase_df['k_off'] = k_off_gyrase
gyrase_df.to_csv('calibration_gyrase.csv', index=False, sep=',')

# And the one with all the values
both_df = pd.DataFrame(columns=['k_on_topoI', 'k_off_topoI', 'k_cat_topoI', 'k_on_gyrase', 'k_off_gyrase', 'k_cat_gyrase', 'sigma0_gyrase'])
both_df['k_cat_topoI'] = best_df['k_cat_topoI']
both_df['k_cat_gyrase'] = best_df['k_cat_gyrase']
both_df['k_on_topoI'] = k_on_topoI
both_df['k_off_topoI'] = k_off_topoI
both_df['k_on_gyrase'] = k_on_gyrase
both_df['k_off_gyrase'] = k_off_gyrase
both_df['sigma0_gyrase'] = sigma_f_gyrase
both_df.to_csv(file_out + '.csv', index=False, sep=',')
