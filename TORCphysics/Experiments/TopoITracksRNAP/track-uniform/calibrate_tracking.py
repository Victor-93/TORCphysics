import numpy as np
import multiprocessing
import pandas as pd
from TORCphysics import parallelization_tools as pt

#TODO:
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
n_simulations = 10#24#8 #96 # 120
tests = 2  # number of tests for parametrization


# Simulation conditions
# --------------------------------------------------------------
dt = 0.25
initial_time = 0
final_time = 100 #500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Reference - It is the reference density of topos when there is no gene that we will use to calculate the
#             fold enrichment.
# --------------------------------------------------------------
reference_topoI_file = '../noRNAP/histogram_topoI.txt'

# Circuit initial conditions
# --------------------------------------------------------------
circuit_filename = '../circuit.csv'
sites_filename = '../sites.csv'
enzymes_filename = '../enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'topoIRNAPtrack'
series = True
continuation = False

initial_sigma = -0.044 # The initial superhelical density
                       # Let's asume it starts at a value where topos are equilibrated, so we assume the steady state.

# Models to calibrate
# -----------------------------------
# Topoisomerase I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'TopoIRecognitionRNAPTracking'

# RANGES FOR RANDOM SEARCH
# -----------------------------------
RNAP_dist_min = 20
RNAP_dist_max = 500
fold_change_min = .1
fold_change_max = 50

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
                         'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': initial_sigma,
                         'DNA_concentration': 0.0}

    # Variation dictionaries
    # ------------------------------------------

    # Topoisomerase I
    name = topoI_name
    object_type = topoI_type
    binding_model_name = topoI_binding_model_name
    binding_oparams = {'k_on': topoI_params['k_on_topoI'], 'width': topoI_params['width_topoI'],
                       'threshold': topoI_params['threshold_topoI'],
                       'RNAP_dist': params['RNAP_dist'],
                       'fold_change': params['fold_change']}

    topoI_variation = {'name': name, 'object_type': object_type,
                       'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams}

    # Create lists of conditions for each system
    # ------------------------------------------

    # Global dictionaries
    global_dict_list = [global_dict]

    # List of lists of variations
    variations_list = [ [topoI_variation] ]

    # Arrays with position densities to calculate fold change
    list_reference = [reference_topoI]

    # TODO: Aqui me quede. Falta terminar de proponer la function, hacer el hyper opt, crear el Tracking
    #  calibration tool, y definir o usar el parallelization tool necesario.
    my_objective = 0 # Borra esto


    # Finally, run objective function. run_objective_function will process our conditions
    # ------------------------------------------
#    my_objective, simulation_superhelicals = tct.run_objective_function(global_dict_list=global_dict_list,
#                                                                        variations_list=variations_list,
#                                                                        exp_superhelicals=list_sigmas,
#                                                                        n_simulations=n_simulations)
    return my_objective

# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------

# Load reference files
reference_topoI = np.loadtxt(reference_topoI_file)

# Load topoI params
topoI_params = pd.read_csv(reference_topoI).to_dict()



