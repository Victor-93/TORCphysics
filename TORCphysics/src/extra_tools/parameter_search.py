import numpy as np
import multiprocessing
import multiprocessing.pool
from TORCphysics import parallelization_tools as pt
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
from TORCphysics import Circuit
from TORCphysics import analysis as ann
import pandas as pd

# The porpuse of this script is to help run calibration processes where we want to find optimum parametrizations
# that reproduce certain behaviours.

# Calibrate according rates given a reference system
def calibrate_w_rate(info_list, target_dict, n_simulations, additional_results=False):

    # Prepare variables
    objective=0.0
    n_systems = len(info_list)
    output_list = []  # Let's return it as a list as well

    # Create a multiprocessing pool
    pool = multiprocessing.Pool()

    # Run simulations
    # --------------------------------------------------------------
    # Go through each system, run simulations in parallel, collect outputs, repeat
    for i in range(n_systems):

        # This contains all the info we need to run the simulation and apply variations
        system = info_list[i]

        # We need a list of items, so the pool can pass each item to the function
        Items = []
        for simulation_number in range(n_simulations):
            g_dict = dict(system['global_conditions'])
            g_dict['n_simulations'] = simulation_number
            Item = {'global_conditions': g_dict, 'variations': system['variations']}
            Items.append(Item)

        # Run in parallel
        # ----------------------------
        # Run simulations in parallel within this subset
        pool_results = pool.map(pt.single_simulation_w_variations_return_dfs, Items)

        # Process transcripts - Serial
        # ----------------------------
        # Collect results (is it better to do it in serial or parallel? What causes more overhead?)
        transcripts = 0
        for result in pool_results:
            sites_df = result['sites_df']
            mask = sites_df['name'] == target_dict['reporter']
            unbinding_event = sites_df[mask]['unbinding'].to_numpy()

            # Calculate number of transcripts produced by the reporter  (no need to do the rate as the time will be canceled out)
            transcripts += np.sum(unbinding_event[:])
        system['transcripts'] = float(transcripts)

        # Here process additional stuff if
        # TODO: Create the processing function and launch calibration!

    # Objective function part
    # --------------------------------------------------------------
    # We need to calculate the relative rate and add to the objective function
    ref_transcript = [d for d in info_list if d['name'] == target_dict['reference_system']][0]['transcripts']

    if ref_transcript <= 0:
        objective += 100  #Something big because it'll give inf or NaN
    else:
        for i in range(n_systems):
            system = info_list[i]
            relative_rate = system['transcripts']/ref_transcript
            system['relative_rate'] = relative_rate
            objective+= (system['reference']-relative_rate)**2

    return objective, output_list

    #if additional_results:
    #    return objective, output_list
    #else:
    #    return objective
