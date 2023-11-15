import numpy as np
import multiprocessing
from TORCphysics import parallelization_tools as pt


# ----------------------------------------------------------------------------------------------------------------------
# Description
# ----------------------------------------------------------------------------------------------------------------------
# This module serves as tools for doing the topoisomearase calibration.

# ----------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------
def Michael_Menten_equation(vmax, KM, S):
    return vmax * S / (KM + S)


def integrate_MM(vmax, KM, substrate0, product0, frames, dt):
    substrate_array = np.zeros(frames)
    product_array = np.zeros(frames)
    substrate_array[0] = substrate0
    product_array[0] = product0

    substrate = substrate0
    product = product0

    for k in range(1, frames):
        v = Michael_Menten_equation(vmax=vmax, KM=KM, S=substrate)
        product = product + v * dt
        substrate = substrate - v * dt
        substrate_array[k] = substrate
        product_array[k] = product
    return substrate_array, product_array


def rescale_product_to_sigma(product, sigma_min, sigma_max):
    #    current_min
    #    frames = len(product)
    #    sigma = np.zeros(frames)
    #    sigma = product-np.min(product)
    #    sigma = sigma/np.max(sigma) # Normalized
    #    d = sigma_f - sigma_i
    #    sigma = (sigma - sigma_i)/1#abs(sigma_f-sigma_i)

    current_min = np.min(product)
    current_max = np.max(product)
    range_ = current_max - current_min
    desired_range = sigma_max - sigma_min
    sigma = ((product - current_min) / range_) * desired_range + sigma_min
    return sigma


# global_dict = Dict with global simulation conditions
# variations_list = A list with a list of variations to implement to the enzymes, environmentals or sites.
# initial_substrates = A list with DNA DNA concentration
# exp_superhelicals = list with values of experimental superhelical densities
# n_simulations = how many simulations to launch.
def run_objective_function(global_dict, variations_list, initial_substrates, exp_superhelicals, n_simulations):
    # Let's run experiments for the substrate concentrations
    my_objective = 0.0
    simulation_superhelicals = []
    for s, substrate0 in enumerate(initial_substrates):
        exp_superhelical = exp_superhelicals[s]  # Experimental superhelical densities
        global_dict['DNA_concentration'] = substrate0  # DNA concentration

        # Let's create an Item to pass the conditions to the simulation
        Item = {'global_conditions': global_dict, 'variations': variations_list}

        # But we actually need a list of items, so the pool can pass each item to the function
        Items = []
        for simulation_number in range(n_simulations):
            g_dict = dict(global_dict)
            g_dict['n_simulations'] = simulation_number
            Item = {'global_conditions': g_dict, 'variations': variations_list}

            Items.append(Item)

        # Create a multiprocessing pool
        pool = multiprocessing.Pool()
        pool_results = pool.map(pt.single_simulation_calibration_w_supercoiling, Items)

        my_supercoiling = np.zeros((global_dict['frames'], n_simulations))
        for i, sigma in enumerate(pool_results):
            my_supercoiling[:, i] = sigma[:-1]
        mea = np.mean(my_supercoiling, axis=1)
        my_objective += np.sum(np.square(np.mean(my_supercoiling, axis=1) - exp_superhelical))
        simulation_superhelicals.append(mea)

    my_objective = my_objective + 0.0
    return my_objective, simulation_superhelicals
