import numpy as np
from TORCphysics import Circuit
import random
import sys
import pickle
import objective_single_molecule_experiment as osme
from TORCphysics import binding_model as bm
from custom_gyrase import GyraseCyclesEffect, GyraseCyclesUnbinding, GyraseCyclesForce
from TORCphysics import Site
import copy
import parallelization_module as pm

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# The idea is that this script re-runs Bustamante's experiment with any model we want and
# any parameterisation.

# We want to run the following scenario: Gyrase acting on relaxed DNA

# And we want to analyse how much twist they add throughout the simualtions...

# Make sure to run the simulations with the same conditions.
# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
output_file = 'gyrase_re-run'
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 0.1 #0.25
initial_time = 0
final_time = 10800
#final_time = 10800 // 2
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
n_sims = 10#100  # This is the total

# For the simulation
# Concentrations in nM
gyrase_concentration = 10.0
#gyrase_concentration = 44.6

# Superhelical values (sigma) for each case
sigma0_RX = 0.0  # Superhelical density of relaxed DNA
# At this value the torque is too strong.

# Models to calibrate to calibrate
# -----------------------------------
gyrase_binding_model = bm.PoissonBinding()
gyrase_effect_model = GyraseCyclesForce()
gyrase_unbinding_model = GyraseCyclesUnbinding()
gyrase_k_cat=20

# You can copy them or load a csv and convert to a dict
gyrase_params = {'k_go': 0.4555795727171365, 'k_off': 0.010492818755306462,
                 'k_on': 0.0909692358724055, 'k_unwrap': 0.487789224391846,
                 'k_wrap': 0.9337534202730736, 'x_wrap': 16.25429900027436,
                 'k_cat': gyrase_k_cat}


single_molecule_size = 2200
Forces = [0.35, 0.6, 0.8, 1.0, 1.3] # pN
output_prefix = 'test0'
series = True
continuation = False

# ----------------------------------------------------------------------------------------------------------------------
# Simulations
# ----------------------------------------------------------------------------------------------------------------------

# The loss and total_rotations will be the outputs
loss = 0
total_rotations = []

# For gyrase Relaxed ------------------------------------
# Let's go through forces
for fi, my_force in enumerate(Forces):  # fi is the index

    # Sort parameters
    gyrase_params['force'] = my_force
    gyrase_params['k_cat'] = gyrase_k_cat

    gyrase_bm = gyrase_binding_model.__class__(**gyrase_params)
    gyrase_em = gyrase_effect_model.__class__(**gyrase_params)
    gyrase_unbm = gyrase_unbinding_model.__class__(**gyrase_params)

    # Define the sites for both enzymes
    # ------------------------------------------
    gyrase_site = Site(site_type='gyrase', name='gyrase_s', start=1000, end=0,
                       binding_model=copy.deepcopy(gyrase_bm))

    # ------------------------------------------
    global_dict_gyrase = {'frames': frames, 'dt': dt, 'superhelicity': sigma0_RX, 'system_name': 'gyrase_RX',
                           'size': single_molecule_size, 'structure': 'linear',
                           'gyrase_concentration': gyrase_concentration,
                           'topoI': False, 'gyrase': True,
                           'gyrase_bm': gyrase_bm, 'gyrase_em': gyrase_em, 'gyrase_unbm': gyrase_unbm,
                           'sites': [gyrase_site]
                           }

    enzymes_df_list = []
    sites_df_list = []
    # Run experiments
    # ----------------------------------------------------
    for simulation_number in range(n_sims):  # Basically, we just pass the info on each global dict,

        global_dict_gyrase['n_sim'] = simulation_number

        # This function runs Bustamantes experiment and applies all the conditions in global_dict_gyrase
        enzymes_df, sites_df, environment_df = pm.single_molecule_simulation_return_dfs(global_dict_gyrase)

        # Append dataframes to the lists.
        enzymes_df_list.append(enzymes_df)
        sites_df_list.append(sites_df)

    error, rotations = osme.Bustamante_objective(
        [{'sites_df_list': sites_df_list, 'enzymes_df_list': enzymes_df_list}],
        my_force)

    loss += error
    total_rotations.append(rotations)

output_dict = {'loss': loss, 'rotations_distribution': total_rotations}
# We can save python objects in binary with pickle.
with open(output_file+'.pkl', 'wb') as file:
    pickle.dump(output_dict, file)
