import sys
#sys.path.append("/mnt/parscratch/users/username")  # For stanage
from custom_gyrase import GyraseCyclesEffect, GyraseCyclesUnbinding
import numpy as np
from TORCphysics import Circuit, TopoisomeraseIEnvironment, GyraseEnvironment
import copy
import multiprocessing  # This one is for parallelizing the process

# This module contains tools for parallelizing the process

# ----------------------------------------------------------------------------------------------------------------------
# We assign the effect and unbinding model.
# For the Binding Model, we use the one we already had, the GyraseRecognition Model (this may change in the future).
def assign_GyraseCycles_models(circuit, **oparams):

    # Binding Model could be assigned in a similar way, but let's not do it yet.

    # Assign effect model
    gyrase_environment = [d for d in circuit.environmental_list if d.name == 'gyrase'][0]
    gyrase_effect_model = GyraseCyclesEffect(**oparams)
    gyrase_environment.effect_model = gyrase_effect_model
    gyrase_environment.effect_model_name = gyrase_effect_model.__class__.__name__
    gyrase_environment.effect_model_oparams = gyrase_effect_model.oparams
    gyrase_environment.effect_model_oparams = None

    # Unbinding Model
    gyrase_environment = [d for d in circuit.environmental_list if d.name == 'gyrase'][0]
    gyrase_unbinding_model = GyraseCyclesUnbinding(**oparams)
    gyrase_environment.unbinding_model = gyrase_unbinding_model
    gyrase_environment.unbinding_model_name = gyrase_unbinding_model.__class__.__name__
    gyrase_environment.unbinding_model_oparams = gyrase_unbinding_model.oparams
    gyrase_environment.unbinding_model_oparams = None

def extract_params(params_dict, name):
    """
    Extract ALL parameters from a dictionary whose key ends with "_<name>".

    Example:
        params = {"k_on_topoI": 1, "k_off_topoI": 2, "Km_topoI": 3}
        extract_params(params, "topoI")
        → {"k_on": 1, "k_off": 2, "Km": 3}
    """
    suffix = f"_{name}"
    extracted = {}

    for key, value in params_dict.items():
        if key.endswith(suffix):
            new_key = key[: -len(suffix)]  # remove the "_name" suffix
            extracted[new_key] = value

    return extracted


# This function will run objective functions for different systems with different conditions.
# The final objective function will be the sum of all objective functions.
# global_dict_list = List of dictionaries with global simulation conditions
# exp_superhelicals = list with arrays of superhelical densities for each system/experiment
# n_simulations = how many simulations to launch per system.
def run_objective_function_exp1(global_dict_list, exp_superhelicals, n_simulations, parallelization=False):
    n_systems = len(global_dict_list)  # number of systems, e.g., topoI alone, gyrase alone, both on rx DNA, both on sc DNA
    # Initialize outputs
    my_objective = 0.0
    simulation_superhelicals = []

    # Iterate over each system
    for n in range(n_systems):

        frames = global_dict_list[n]['frames']
        sname = global_dict_list[n]['system_name']

            # We need a list of items, so the pool can pass each item to the function - right now is single cpu but its ok
        Items = []
        for simulation_number in range(n_simulations): # Basically, we just pass the info on each global dict,
                                                       # but assign a simulation number
            g_dict = copy.deepcopy(global_dict_list[n])
            g_dict['n_sim'] = simulation_number
            Items.append(g_dict)

        my_supercoiling = np.zeros((frames, n_simulations)) # Initialize array to process superhelical density

        if parallelization:
            # Run process in parallel
            with multiprocessing.Pool() as pool:
                pool_results = pool.map(single_simulation_topo_calibration_r_supercoiling, Items)

            # Process superhelical densities to calculate objective function
            my_supercoiling = np.zeros((frames, n_simulations)) # Initialize array to process superhelical density
            for i, sigma in enumerate(pool_results):
                my_supercoiling[:, i] = sigma[:-1]

        else:
            for i in range(n_simulations):  # Run multiple simulations
                my_supercoiling[:, i] = single_simulation_topo_calibration_r_supercoiling(Items[i])[:-1]

        mean = np.mean(my_supercoiling, axis=1)
        std = np.std(my_supercoiling, axis=1)
        current_objective = np.sum(np.square(np.mean(my_supercoiling, axis=1) - exp_superhelicals[n]))

        # Save average superhelical densities
        simulation_superhelicals.append({'name': sname, 'mean':mean, 'std': std})

        # The objective (or error) is the sum of errors for each system
        my_objective = my_objective + current_objective

    return my_objective, simulation_superhelicals

# Similar to the run_objective_function, but it runs a more general simulation of topoisomerases
# with specific binding sites. It currently doesnt have a processing step.
def run_experiment_2(global_dict_list, n_simulations, parallelization=False):
    n_systems = len(global_dict_list)  # number of systems, e.g., topoI alone, gyrase alone, both on rx DNA, both on sc DNA
    # Initialize outputs
    output_list = []  # Here we will collect the output.

    # Iterate over each system
    for n in range(n_systems):

        # Extract info we need
        sname = global_dict_list[n]['system_name']

            # We need a list of items, so the pool can pass each item to the function - right now is single cpu but its ok
        Items = []
        for simulation_number in range(n_simulations): # Basically, we just pass the info on each global dict,
                                                       # but assign a simulation number
            g_dict = copy.deepcopy(global_dict_list[n])
            g_dict['n_sim'] = simulation_number
            Items.append(g_dict)

        # Prepare output lists
        enzymes_df_list = []
        sites_df_list = []

        # RUN
        if parallelization:
            # Run process in parallel
            with multiprocessing.Pool() as pool:
                pool_results = pool.map(single_molecule_simulation_return_dfs, Items)

            # IF we want to do some pre-processing or collect other stuff, we would probably do it here.

            # Process each pool result
            for i, pool in enumerate(pool_results):
                enzymes_df_list.append(pool[0])
                sites_df_list.append(pool[1])
        else:
            for i in range(n_simulations):  # Run multiple simulations
                enzymes_df, sites_df, environment_df = single_molecule_simulation_return_dfs(Items[i])

                # IF we want to do some pre-processing or collect other stuff, we would probably do it here.

                # Collect data we need
                enzymes_df_list.append(enzymes_df)
                sites_df_list.append(sites_df)

        # Save average superhelical densities
        output_list.append({'name':sname, 'enzymes_df_list': enzymes_df_list, 'sites_df_list': sites_df_list})

    return output_list

def single_simulation_topo_calibration_r_supercoiling(gd):
    """
    Function to run a simulation of topoisomerases.  It only returns the global superhelical
    densities which are used to calculate the objective function (error).

    Inputs
     gd: a global dictionary with simulation conditions.

    Returns
     supercoiling: global superhelical density
    """
    # gd = global dictionary with all information we need

    # Initialize circuit with topoisomerases
    # --------------------------------------------------------------------------
    my_circuit = Circuit(name='topos_'+str(gd['n_sim']), structure=gd['structure'],
                           superhelicity=gd['superhelicity'], size=gd['size'],
                           frames=gd['frames'], dt=gd['dt'])
    # Add topoisomerases
    if gd['topoI']:
        TopoisomeraseIEnvironment(binding_model=copy.deepcopy(gd['topoI_bm']),
                                  effect_model = copy.deepcopy(gd['topoI_em']),
                                  unbinding_model = copy.deepcopy(gd['topoI_unbm']),
                                  concentration=gd['topoI_concentration'],
                                  circuit=my_circuit)
    if gd['gyrase']:
        GyraseEnvironment(binding_model=copy.deepcopy(gd['gyrase_bm']),
                          effect_model=copy.deepcopy(gd['gyrase_em']),
                          unbinding_model=copy.deepcopy(gd['gyrase_unbm']),
                          circuit=my_circuit)

    # Slightly change the seed
    my_circuit.seed = my_circuit.seed + gd['n_sim']
    my_circuit.rng = np.random.default_rng(my_circuit.seed)

    # Run simulation but only collecting global superhelical densities
    supercoiling = my_circuit.run_return_global_supercoiling()
    return supercoiling

def single_molecule_simulation_return_dfs(gd):
    """
    Runs a simulation for the single_molecule experiment of topoisomerases.
    Simulation conditions are provided by the global dict gd.
    Here, we assume topoisomerases bind to specific sites, not everywhere.

    Inputs
     gd: global dictionary

    Returns
     enzymes_df: enzymes dataframe
     sites_df:  sites dataframe
     environment_df:  environment dataframe
    """
    # Initialize circuit with topoisomerases
    # --------------------------------------------------------------------------
    my_circuit = Circuit(name='topos_'+str(gd['n_sim']), structure=gd['structure'],
                           superhelicity=gd['superhelicity'], size=gd['size'],
                           frames=gd['frames'], dt=gd['dt'])

    # Add sites to the circuit
    for site in gd['sites']:
        my_circuit.add_custom_Site(copy.deepcopy(site))

    # Add topoisomerases
    if gd['topoI']:
        TopoisomeraseIEnvironment(binding_model=copy.deepcopy(gd['topoI_bm']), binding_model_name=None,
                                  effect_model = copy.deepcopy(gd['topoI_em']),
                                  unbinding_model = copy.deepcopy(gd['topoI_unbm']),
                                  concentration=gd['topoI_concentration'],
                                  # site_list=utils.site_match_by_name(my_circuit.site_list, 'topoI'),
                                  site_type='topoI',
                                  circuit=my_circuit
                                  )
    if gd['gyrase']:
        GyraseEnvironment(binding_model=copy.deepcopy(gd['gyrase_bm']), binding_model_name=None,
                          effect_model=copy.deepcopy(gd['gyrase_em']),
                          unbinding_model=copy.deepcopy(gd['gyrase_unbm']),
                          # site_list=utils.site_match_by_name(my_circuit.site_list, 'gyrase'),
                          site_type='gyrase',
                          circuit=my_circuit)

    # Slightly change the seed
    my_circuit.seed = my_circuit.seed + gd['n_sim']
    my_circuit.rng = np.random.default_rng(my_circuit.seed)

    # Run simulation
    enzymes_df, sites_df, environment_df = my_circuit.run_return_dfs()
    return enzymes_df, sites_df, environment_df

if __name__ == "__main__":

    print('Hi')