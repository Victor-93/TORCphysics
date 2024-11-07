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
        additional_dict = {}

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

        # TODO: Let's add prod_rates, cross_correlations, local/global supercoiling levels  and  KDEs?
        # Process transcripts - Serial
        # ----------------------------
        # Collect results (is it better to do it in serial or parallel? What causes more overhead?)
        transcript_list = []
        for result in pool_results:
            environmental_df = result['environmental_df']
            mask = environmental_df['name'] == target_dict['reporter']

            if len(environmental_df[mask]['concentration']) == 0:
                transcript = 0
            else:
                transcript = environmental_df[mask]['concentration'].iloc[-1]
            transcript_list.append(transcript)

            # sites_df = result['sites_df']
            # mask = sites_df['name'] == target_dict['reporter']
            # unbinding_event = sites_df[mask]['unbinding'].to_numpy()

            # Calculate number of transcripts produced by the reporter  (no need to do the rate as the time will be canceled out)
            # transcripts += np.sum(unbinding_event[:])

        # Convert to a NumPy arraypool_results
        transcripts_array = np.array(transcript_list)

        # Calculate mean and standard deviation
        mean_transcripts = np.mean(transcripts_array)
        std_transcripts = np.std(transcripts_array)

        transcripts = np.array([mean_transcripts, std_transcripts])

        system['transcripts'] = transcripts

        # Here process additional stuff if
        # TODO: Create the processing function and launch calibration! - This would be parameters related to the dfs
        #  Maybe think about this when you actually need it, because you need to pass additional information like x_gene
        #  etc... Maybe build those arrays some part in the beggining of this function if additional_results=True
        if additional_results:   # Note I add it here to start giving shape to the code, in case in the future
                                 # I want to add the KDEs or the processing bit, so it can

            # Total transcripts
            additional_dict['transcripts'] = transcripts
            output_list.append(additional_dict)

            # System info
            output_list[i]['name'] = system['name']
            output_list[i]['despcription'] = system['description']
            output_list[i]['bacterium'] = system['bacterium']
            output_list[i]['promoter'] = system['promoter']
            output_list[i]['strain'] = system['strain']

            output_list[i]['reference'] = np.array([system['reference'], system['reference_std']])

            # Get Supercoiling and related df quantities
            # TODO: AQUIMEQUEDE: You have somehow to save local sigma, global sigma per system
            # Maybe the positions as well
            # --------------------------------------------------------------------
            for result in pool_results:
                sites_df = result['sites_df']
                enzymes_df = result['enzymes_df']
                site_mask = sites_df['name'] == target_dict['reporter']
                mask_circuit = sites_df['type'] == 'circuit'

                # Collect measurements
                local_df = sites_df[site_mask]#['superhelical']
                local_superhelical = local_df['superhelical'].to_numpy()
                # local_superhelical = sites_df[site_mask]['superhelical'].to_numpy()
                global_superhelical = sites_df[mask_circuit]['superhelical'].to_numpy()



    # Objective function part
    # --------------------------------------------------------------
    # We need to calculate the relative rate and add to the objective function
    ref_transcript = [d for d in info_list if d['name'] == target_dict['reference_system']][0]['transcripts'][0]

    if ref_transcript <= 0:
        objective += 100  #Something big because it'll give inf or NaN
        if additional_results:
            output_list[i]['objective'] = 0
            output_list[i]['relative_rate'] = 0
    else:
        for i in range(n_systems):
            system = info_list[i]
            relative_rate = system['transcripts']/ref_transcript
            system['relative_rate'] = relative_rate
            system_objective = (system['reference']-relative_rate[0])**2
            objective+= system_objective

            # Add values to additional outputs
            if additional_results:
                output_list[i]['objective'] = system_objective
                output_list[i]['relative_rate'] = relative_rate

    return objective, output_list



    #if additional_results:
    #    return objective, output_list
    #else:
    #    return objective

#def process_additional_results(item):
#    # Global measurements

#    mask_circuit = sites_df['type'] == 'circuit'
#    global_superhelical = sites_df[mask_circuit]['superhelical'].to_numpy()

#    additional_dict['global_superhelical'] = global_superhelical
