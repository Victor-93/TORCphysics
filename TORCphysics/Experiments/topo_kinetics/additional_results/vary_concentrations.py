import numpy as np
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
import pickle

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# We launch simulations in parallel varying concentrations. We will output three arrays: Global superhelical level
# from simulations, steady state (#enzymes), and the global superhelical from the MM equation.
# We will output these arrays in the form of [x,y,ys] where x is time, y the quantity of interest, and ys the
# standard deviation. Notice that the superhelical density from the MM equation does not have a ys.

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)

# For parallelization and calibration
n_simulations = 24#200#100

# Concentrations in nM
DNA_concentration = 0.75
gyrase_concentration_min = 0.0
gyrase_concentration_max = 100.0
gyrase_concentration_step = 2.0
topoI_concentration_min = 0.0
topoI_concentration_max = 100.0
topoI_concentration_step = 2.0

dt = 1.0 #0.25
initial_time = 0
final_time = 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'concentrations-01-dt'+str(dt)

# For the simulation
circuit_filename = 'circuit.csv'
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
environment_filename = 'environment_small.csv'

output_prefix = 'test0'
series = True
continuation = False

recognition_path = '../'
recognition_params_file = 'avg_dt'+str(dt)+'.csv'

# Superhelical values (sigma) for each case
sigma_0_topo = -0.11#0.075  # Approximately -20 supercoils according the paper
sigma_0_gyrase = 0.0  # We suppose this one.
sigma_f_gyrase = -0.11  # We also assume this one, which is the maximum at which gyrase acts.
# At this value the torque is too strong.


# Models to use
# -----------------------------------
# Topoisomerase I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'TopoIRecognition'
topoI_effect_model_name = 'TopoILinear'
topoI_unbinding_model_name = 'PoissonUnBinding'
topoI_params = '../topoI_rec_avg_dt'+str(dt)+'.csv'

# Gyrase
gyrase_name = 'gyrase'
gyrase_type = 'environmental'
gyrase_binding_model_name = 'GyraseRecognition'
gyrase_effect_model_name = 'GyraseLinear'
gyrase_unbinding_model_name = 'PoissonUnBinding'
gyrase_params = '../gyrase_rec_avg_dt'+str(dt)+'.csv'

# MM kinetics
K_M_topoI = 1.5
k_cat_topoI = .0023
K_M_gyrase = 2.7
k_cat_gyrase = .0011

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# MM functions
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------
# Build experimental curve for TOPO I
# -----------------------------------------------------
# Kinetics: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI
# Product = Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs; which initially is the same as the DNA concentration
# Integrate MM kinetics
def topoI_Rx(topoI_c):

    # MM kinetics
    v_max_topoI = k_cat_topoI * topoI_c

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
    return sigma

# -----------------------------------------------------
# Build experimental curve for Gyrase
# -----------------------------------------------------
# Kinetics: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Product = Supercoiled DNA
# Substrate = Relaxed DNA; which initially is the same as the DNA concentration
def gyrase_Sc(gyrase_c):

    # MM kinetics
    v_max_gyrase = k_cat_gyrase * gyrase_c

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
    return sigma

# -----------------------------------------------------
# -----------------------------------------------------
# Build experimental curve for system with both Topo I and Gyrase (from Sc state)
# -----------------------------------------------------
# Kinetics Gyrase: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Kinetics Topoisomerase: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI
def both_Sc(topoI_c, gyrase_c):

    # MM kinetics
    v_max_topoI = k_cat_topoI * topoI_c
    v_max_gyrase = k_cat_gyrase * gyrase_c

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
    return sigma

# -----------------------------------------------------
# Build experimental curve for system with both Topo I and Gyrase (from Rx state)
# -----------------------------------------------------
def both_Rx(topoI_c, gyrase_c):

    # MM kinetics
    v_max_topoI = k_cat_topoI * topoI_c
    v_max_gyrase = k_cat_gyrase * gyrase_c

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
    return sigma
# ----------------------------------------------------------------------------------------------------------------------

# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the mean superhelical
# density for each substrate concentration

def objective_function(params):
    gyrase_concentration = params['gyrase_concentration']
    topoI_concentration = params['topoI_concentration']

    # We need to prepare the inputs.
    # This time we have two different systems:
    # Both enzymes acting on supercoiled/Relaxed DNA

    # Global dictionaries
    # ------------------------------------------
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
    if 'Poisson' in binding_model_name:
        binding_oparams = {'k_on': float(params['k_on_topoI'][0])}
    else:
        binding_oparams = {'k_on': float(params['k_on_topoI'][0]), 'width': float(params['width_topoI'][0]),
                           'threshold': float(params['threshold_topoI'][0])}
    effect_model_name = topoI_effect_model_name
    effect_oparams = {'k_cat': float(params['k_cat_topoI'][0])}
    unbinding_model_name = topoI_unbinding_model_name
    unbinding_oparams = {'k_off': float(params['k_off_topoI'][0])}
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
    if 'Poisson' in binding_model_name:
        binding_oparams = {'k_on': float(params['k_on_gyrase'][0])}
    else:
        binding_oparams = {'k_on': float(params['k_on_gyrase'][0]), 'width': float(params['width_gyrase'][0]),
                           'threshold': float(params['threshold_gyrase'][0])}
    effect_model_name = gyrase_effect_model_name
    effect_oparams = {'k_cat': float(params['k_cat_gyrase'][0]), 'sigma0': float(params['sigma0_gyrase'][0])}
    unbinding_model_name = gyrase_unbinding_model_name
    unbinding_oparams = {'k_off': float(params['k_off_gyrase'][0])}
    concentration = gyrase_concentration  # / mol_concentration  # Because this is the reference.

    gyrase_variation = {'name': name, 'object_type': object_type,
                        'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                        'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                        'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                        'concentration': concentration}

    # Create lists of conditions for each system
    # ------------------------------------------

    # Global dictionaries
    global_dict_list = [global_dict_both_sc, global_dict_both_rx]

    # List of lists of variations
    variations_list = [
                       [topoI_variation, gyrase_variation],
                       [topoI_variation, gyrase_variation]
                       ]

    # Arrays with global superhelical densities
    both_sigma_sc = both_Sc(topoI_concentration, gyrase_concentration)
    both_sigma_rx = both_Rx(topoI_concentration, gyrase_concentration)
    list_sigmas = [both_sigma_sc, both_sigma_rx]

    # Finally, run objective function. run_objective_function will process our conditions
    # ------------------------------------------
    my_objective, output_dict = tct.test_steady_topos(global_dict_list=global_dict_list,
                                                                        variations_list=variations_list,
                                                                        exp_superhelicals=list_sigmas,
                                                                        n_simulations=n_simulations)
    return my_objective, output_dict, list_sigmas

# ----------------------------------------------------------------------------------------------------------------------
# Process objective function - but it is not really an objective function in this case because we are not
#                              calibrating anything
# ----------------------------------------------------------------------------------------------------------------------
# Create concentrations arrays
gyrase_vars = np.arange(gyrase_concentration_min, gyrase_concentration_max + gyrase_concentration_step, gyrase_concentration_step)
topoI_vars = np.arange(topoI_concentration_min, topoI_concentration_max + topoI_concentration_step, topoI_concentration_step)

total_cases = len(gyrase_vars) * len(topoI_vars)
print("Running " + str(n_simulations) + " simulations.")
print("for a total of " + str(total_cases) + " cases.")
output_list = []  #List with outputs
for gyrase_i, gyrase_concentration in enumerate(gyrase_vars):
    for topoI_i, topoI_concentration in enumerate(topoI_vars):

        # Set up input
        params_file = recognition_path+recognition_params_file
        params_dict = pd.read_csv(params_file).to_dict()
        params_dict['gyrase_concentration'] = gyrase_concentration
        params_dict['topoI_concentration'] = topoI_concentration

        # Run the objective function
        objective, sim_output, sigmas = objective_function(params=params_dict)

        # Prepare output dict result_dict
        result_dict = {}
        result_dict['gyrase_concentration'] = gyrase_concentration
        result_dict['topoI_concentration'] = topoI_concentration

        # Collect results
        result_dict['objective'] = objective

        # Because we are running experiments from a Rx and Sc state, we separate each run into these two outputs
        output_sc = sim_output[0]
        output_sc['MM_supercoiling'] = sigmas[0]
        output_rx = sim_output[1]
        output_rx['MM_supercoiling'] = sigmas[1]
        result_dict['both_Sc'] = output_sc
        result_dict['both_Rx'] = output_rx

        output_list.append(result_dict)

# Save the dictionary to a file
with open(file_out + '.pkl', 'wb') as file:
    pickle.dump(output_list, file)

