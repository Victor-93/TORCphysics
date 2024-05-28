import numpy as np
import multiprocessing
from TORCphysics import parallelization_tools as pt
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
from TORCphysics import Circuit
import pandas as pd


# ----------------------------------------------------------------------------------------------------------------------
# Description
# ----------------------------------------------------------------------------------------------------------------------
# This module serves as tools for doing the topoisomearase calibration.

# Params
buf_size = 10  # Used in calculating correlation. Is the number of data points ignored at the ends


# ----------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------
def Michael_Menten_equation(vmax, KM, S):
    return vmax * S / (KM + S)


# Kinetics: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI
# In this function, topoI acts on supercoiled DNA and produced relaxed DNA
# Substrate = supercoiled_DNA
# product = relaxed DNA
# Supercoiled_0 : Initial concentration of supercoiled DNA
# Relaxed_0 : Initial concentration of relaxed DNA
def integrate_MM_topoI(vmax, KM, Supercoiled_0, Relaxed_0, frames, dt):
    Supercoiled = np.zeros(frames)
    Relaxed = np.zeros(frames)
    Supercoiled[0] = Supercoiled_0
    Relaxed[0] = Relaxed_0

    SC = Supercoiled_0
    RE = Relaxed_0

    for k in range(1, frames):
        v = Michael_Menten_equation(vmax=vmax, KM=KM, S=SC)
        RE = RE + v * dt
        SC = SC - v * dt
        Supercoiled[k] = SC
        Relaxed[k] = RE
    return Supercoiled, Relaxed


# Relaxed: Array with concentration of Relaxed DNA as a function of time.
# DNA_concentration: Parameter with concentration of DNA. Should be equivalent to the initial amount of
# concentration of supercoiled DNA.
# sigma0: Initial level of superhelical density
def topoI_to_sigma(Relaxed, DNA_concentration, sigma0):
    sigma = np.zeros_like(Relaxed)
    n = len(Relaxed)
    for i in range(n):
        sigma[i] = sigma0 - Relaxed[i] * sigma0 / DNA_concentration
    return sigma


# Kinetics: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# In this function, gyrase acts on Relaxed DNA and produces supercoiled DNA
# Product = Supercoiled DNA
# Substrate = Relaxed DNA; which initially is the same as the DNA concentration
# Supercoiled_0 : Initial concentration of supercoiled DNA
# Relaxed_0 : Initial concentration of relaxed DNA
def integrate_MM_gyrase(vmax, KM, Supercoiled_0, Relaxed_0, frames, dt):
    Supercoiled = np.zeros(frames)
    Relaxed = np.zeros(frames)
    Supercoiled[0] = Supercoiled_0
    Relaxed[0] = Relaxed_0

    SC = Supercoiled_0
    RE = Relaxed_0

    for k in range(1, frames):
        v = Michael_Menten_equation(vmax=vmax, KM=KM, S=RE)
        SC = SC + v * dt
        RE = RE - v * dt
        Supercoiled[k] = SC
        Relaxed[k] = RE
    return Supercoiled, Relaxed


# Relaxed: Array with concentration of Relaxed DNA as a function of time.
# DNA_concentration: Parameter with concentration of DNA. Should be equivalent to the initial amount of
# concentration of Relaxed DNA.
# sigma0: Initial level of superhelical density
# sigmaf: Final level of superhelical density
def gyrase_to_sigma(Relaxed, DNA_concentration, sigma0, sigmaf):
    sigma = np.zeros_like(Relaxed)
    n = len(Relaxed)
    for i in range(n):
        sigma[i] = sigmaf + Relaxed[i] * (sigma0 - sigmaf) / DNA_concentration
    return sigma


# Kinetics Gyrase: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Kinetics Topoisomerase: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI
# In this function, both topo I and gyrase are active.
# Gyrase acts on Relaxed DNA and produces supercoiled DNA, while topo I acts on supercoiled DNA and produces relaxed
# DNA
# Supercoiled_0 : Initial concentration of supercoiled DNA
# Relaxed_0 : Initial concentration of relaxed DNA
def integrate_MM_both_T_G(vmax_topoI, vmax_gyrase, KM_topoI, KM_gyrase, Supercoiled_0, Relaxed_0, frames, dt):
    Supercoiled = np.zeros(frames)
    Relaxed = np.zeros(frames)
    Supercoiled[0] = Supercoiled_0
    Relaxed[0] = Relaxed_0

    SC = Supercoiled_0
    RE = Relaxed_0

    for k in range(1, frames):
        v_gyrase = Michael_Menten_equation(vmax=vmax_gyrase, KM=KM_gyrase, S=RE)
        v_topoI = Michael_Menten_equation(vmax=vmax_topoI, KM=KM_topoI, S=SC)
        SC = SC + (v_gyrase - v_topoI) * dt
        RE = RE - (v_gyrase - v_topoI) * dt
        Supercoiled[k] = SC
        Relaxed[k] = RE
    return Supercoiled, Relaxed


# This function only works for when intial  [Relaxed] ~ 0, as I deduced the equation from this condition.
# Relaxed: Array with concentration of Relaxed DNA as a function of time.
# Relaxed_final : Parameter with final concentration of relaxed DNAs.
# DNA_concentration: Parameter with concentration of DNA. Should be equivalent to the initial amount of
# concentration of Relaxed DNA.
# sigma0: Initial level of superhelical density
# sigmaf: Final level of superhelical density
def both_T_G_to_sigma(Relaxed, Relaxed_final, sigma0, sigmaf):
    sigma = np.zeros_like(Relaxed)
    n = len(Relaxed)
    for i in range(n):
        sigma[i] = sigma0 + Relaxed[i] * (sigmaf - sigma0) / Relaxed_final
    return sigma


# This function simply translates the superhelical density from the concentration of
# relaxed DNA curve.
def sigma_to_relaxed(Relaxed, sigmaf, DNA_concentration):
    sigma = np.zeros_like(Relaxed)
    n = len(Relaxed)
    for i in range(n):
        sigma[i] = sigmaf - sigmaf * Relaxed[i] / DNA_concentration
    return sigma


# This function will run objective functions for different systems with different conditions.
# The final objective function will be the sum of all objective functions.
# global_dict_list = List of dictionaries with global simulation conditions
# variations_list = A list with a list of variations to implement to the enzymes, environmentals or sites.
# exp_superhelicals = list with arrays of superhelical densities for each system/experiment
# n_simulations = how many simulations to launch per system.
def run_objective_function(global_dict_list, variations_list, exp_superhelicals, n_simulations):
    n_systems = len(global_dict_list)  # number of systems.

    # Let's run experiments for the substrate concentrations
    my_objective = 0.0
    simulation_superhelicals = []

    for n in range(n_systems):

        # We need a list of items, so the pool can pass each item to the function
        Items = []
        for simulation_number in range(n_simulations):
            g_dict = dict(global_dict_list[n])
            g_dict['n_simulations'] = simulation_number
            Item = {'global_conditions': g_dict, 'variations': variations_list[n]}
            Items.append(Item)

        # Create a multiprocessing pool
        pool = multiprocessing.Pool()
        pool_results = pool.map(pt.single_simulation_calibration_w_supercoiling, Items)

        # Process superhelical densities to calculate objective function
        my_supercoiling = np.zeros((g_dict['frames'], n_simulations))
        for i, sigma in enumerate(pool_results):
            my_supercoiling[:, i] = sigma[:-1]

        mea = np.mean(my_supercoiling, axis=1)
        current_objective = np.sum(np.square(np.mean(my_supercoiling, axis=1) - exp_superhelicals[n]))

        # Save average superhelical densities
        simulation_superhelicals.append(mea)

        #print('system', n, 'simulation', simulation_number, 'my_objective', my_objective)
        my_objective = my_objective + current_objective
    return my_objective, simulation_superhelicals


# Runs the objective function for the topoisomerase I RNAP Tracking model for different system conditions.
# The objective function is calculated as a combination of fold change and correlation between densities of
# topo I and RNAP
# g_dict = global_dict_list = List of dictionaries with global simulation conditions
# variations_list = A list with a list of variations to implement to the enzymes, environmentals or sites.
# reference_arrays = A list of arrays that contain the reference density (position) of topo I in a system
#                    in which there is no gene expression.
# n_simulations = how many simulations to launch per system.
# target_FE = It is the target fold-enrichment we want to achieve for topo I.
# target_CO = Target correlation
# Returns: my_objective (objective value), output_dict - An output in the form of dictionary.
#  This output contains {FE (fold enrichment), correlation (correlation between topoI and RNAP), results}
#   results is a dict, that contains the environmentals kde on the whole system, and on the target gene.
#   It also contains the fold_enrichment curves.
def single_case_RNAPTracking_calibration(global_dict_list, variations_list, reference_dict, n_simulations,
                                         target_dict):
    #n_systems = len(global_dict_list)  # number of systems.

    g_dict = dict(global_dict_list[0])  # Just one system
    # densities_outputs = {}  # This one will contain the outputs in the form of dict, where each enzyme_name is the key
    # that will contain another dict with the results: positions, kde_x, kde_y

    # Names to filter
    enzymes_names = target_dict['enzymes_names']

    # These three dictionaries will be part of our output dict
    kde_gene = {}
    kde_system = {}
    FE_dict = {}
    hists_dict = {}
    for names in enzymes_names:
        kde_gene[names] = None
        kde_system[names] = None
        FE_dict[names] = None
        hists_dict[names] = None
        # densities_outputs[enzymes_names] = None

    #for n in range(n_systems):

    # We need a list of items, so the pool can pass each item to the function
    Items = []
    for simulation_number in range(n_simulations):
        #g_dict = dict(global_dict_list)
        g_dict['n_simulations'] = simulation_number
        Item = {'global_conditions': g_dict, 'variations': variations_list[0]}  #because is list 0
        Items.append(Item)

    # Create a multiprocessing pool
    pool = multiprocessing.Pool()
    pool_results = pool.map(pt.single_simulation_w_variations_return_dfs, Items)

    # Process simulations outputs
    # --------------------------------------------------------------
    # Let's load the circuit, so we can extract the information that we need in an automatic way
    my_circuit = load_circuit(g_dict)

    # Extract Global superhelical density
    # ------------------------------------------------------------------------------
    # Extract info from dataframes
    for j, out_dict in enumerate(pool_results):

        sites_df = out_dict['sites_df']
        mask = sites_df['type'] == 'circuit'
        superhelical = sites_df[mask]['superhelical'].to_numpy()

        simulation_df = pd.DataFrame({'simulation_' + str(j): superhelical})

        if j == 0:
            superhelical_output_df = simulation_df.copy()
        else:
            superhelical_output_df = pd.concat([superhelical_output_df, simulation_df], axis=1).reset_index(drop=True)

    superhelical_mean = superhelical_output_df.mean(axis=1).to_numpy()
    superhelical_std = superhelical_output_df.std(axis=1).to_numpy()

    superhelical_dict = {'mean': superhelical_mean, 'std': superhelical_std}

    # Extract positions
    # ------------------------------------------------------------------------------
    # Get target site
    target_gene = [site for site in my_circuit.site_list if site.name == target_dict['target_gene']][0]
    RNAP_env = [environment for environment in my_circuit.environmental_list if environment.name == 'RNAP'][0]

    # Define x-axes
    x_system = get_interpolated_x(1, my_circuit.size)
    # x_gene = get_interpolated_x(target_gene.start, target_gene.end)
    x_gene = get_interpolated_x(target_gene.start - RNAP_env.size, target_gene.end)

    # And filter results starting by name
    for i, name in enumerate(enzymes_names):

        all_positions = np.empty([1])

        # Extract info from dataframes
        for j, out_dict in enumerate(pool_results):
            enzymes_df = out_dict['enzymes_df']

            # Filter superhelical density
            mask = enzymes_df['name'] == name

            # Position
            # ------------------------------------------------------------------------------
            position = enzymes_df[mask]['position'].to_numpy()

            # position = position[~np.isnan(position)]  # Just in case to remove nans

            #if j == 0:
            #    all_positions = np.copy(position)
            #else:
            #    all_positions = np.concatenate([position, all_positions])
            all_positions = np.concatenate([position, all_positions])

            #s=0+len(position)
        # number of bins for histogram - which we will not plot
        nbins = calculate_number_nbins(my_circuit, name)

        # Calculate histogram
        counts, bin_edges = calculate_histogram(data=all_positions, nbins=nbins)
        hists_dict[name] = {'counts': counts, 'bin_edges': bin_edges}

        # print(name, len(all_positions))
        # Calculate KDE
        kde_x, kde_y = calculate_KDE(data=all_positions, nbins=nbins, scaled=True)

        # Let's interpolate for the gene - For every enzyme, let's interpolate signal in gene region for comparison
        kde_gene[name] = get_interpolated_kde(kde_x, kde_y, x_gene)

        # kde for the system (ony for topos)
        if name != 'RNAP':
            kde_system[name] = get_interpolated_kde(kde_x, kde_y, x_system)

            # Compute the fold-enrichment since we are here (no FE for RNAP)
            FE_dict[name] = kde_system[name] / reference_dict[0][name]

        #print(name, len(all_positions))
        # Calculate the KDE
        # kde = gaussian_kde(all_positions)
        # k#de_x = np.linspace(min(all_positions), max(all_positions), nbins[i])
        # kde_y = kde(kde_x)

        # And save them to our output dict
        # densities_outputs[name] = {'positions': all_positions, 'kde_x': kde_x, 'kde_y': kde_y}

    pool.close()

    # Compute correlation between topoI and RNAP signals along the gene
    # --------------------------------------------------------------------
    buf_size = 10  # Used in calculating correlation. Is the number of data points ignored at the ends
    correlation_matrix = np.corrcoef(kde_gene['topoI'][buf_size:-buf_size], kde_gene['RNAP'][buf_size:-buf_size])
    correlation = correlation_matrix[0, 1]

    # Retrieve the target correlation
    target_CO = target_dict['target_CO']

    # Compute average fold-enrichment FE
    # --------------------------------------------------------------------
    avg_FE = np.mean(FE_dict['topoI'])  # TODO: Add ranges [a:b]? Or not?
    target_FE = target_dict['target_FE']

    # Normalize FE so we can compare with the correlation
    #FE_curve_norm = (FE_dict['topoI'] - np.min(FE_dict['topoI'])) / (np.max(FE_dict['topoI']) - np.min(FE_dict['topoI']))

    # The average of this new normalized curve
    #avg_FE_norm = np.mean(FE_curve_norm)  # If you add ranges up, then here as well

    # And normalize the target_FE in the same way, so we can compare them
    #target_FE_norm = ((target_dict['target_FE'] - np.min(FE_dict['topoI'])) /
    #                  (np.max(FE_dict['topoI']) - np.min(FE_dict['topoI'])))
    #minn = np.min(FE_dict['topoI'])
    #maxx = np.max(FE_dict['topoI'])

    # Note that we did not normalise the correlation as it is already normalised

    # Objective function
    #--------------------------------------------------------------------
    #my_objective = (target_FE_norm - avg_FE_norm) ** 2 + (target_CO - correlation) ** 2
    # Maybe let's not compute the normalized FE as both correlation and FE are comparable in magnitude
    my_objective = (target_FE - avg_FE) ** 2 + (target_CO - correlation) ** 2

    # And prepare output dict
    #--------------------------------------------------------------------
    results_dict = {'kde_gene': kde_gene, 'kde_system': kde_system, 'FE_dict': FE_dict, 'hists_dict': hists_dict,
                    'superhelical_dict': superhelical_dict}
    output_dict = {'FE': avg_FE, 'correlation': correlation, 'results': results_dict}

    return my_objective, output_dict


# This function is similar to single_case_RNAPTracking_calibration()
# The main difference is that it runs a number of M simulations in parallel for N sets (or batches).
# The idea is that for each of these sets, we obtain histograms, kdes, global superhelical densities, FEs and CO.
# KDEs are used to calculate a smooth KDE as well as histograms. From these overall KDEs the correlation between
# topoI and RNAP are calculated, while the FE is calculated as the average of FEs across sets.
# This because the FE do not vary greatly between sets of simulations, while the correlation does. This is the main
# reason of why we calculate many KDEs to get the overalls.
# FE = Fold enrichment
# CO = Correlations
# The results output include: overall histogram, KDEs, FE curves, and global superhelicals with their respective STDs.
# It also includes averaged FEs and correlations, plus the Overall Correlation calculated from overall KDEs.
# The objective function is minimized according the averaged FE of topo I, and the overall correlation.
# global_dict_list should include the number of sets.
def single_case_RNAPTracking_calibration_nsets(global_dict_list, variations_list, reference_dict, n_simulations,
                                               target_dict):
    #n_systems = len(global_dict_list)  # number of systems.

    g_dict = dict(global_dict_list[0])  # Just one system

    n_sets = g_dict['n_sets']  # Extracts number of sets: n_sets

    # Names to filter
    enzymes_names = target_dict['enzymes_names']

    # The these dicts include another dict with mean and STD columns - Everything is dictionary this time.
    kde_gene = {}
    kde_system = {}
    FE_curve = {}
    hists_dict = {}
    FE_val = {}
    correlation_sets = {'mean': None, 'std': None}
    global_supercoiling = {'mean': None, 'std': None}

    for names in enzymes_names:
        kde_gene[names] = {'mean': None, 'std': None}
        kde_system[names] = {'mean': None, 'std': None}
        FE_curve[names] = {'mean': None, 'std': None}
        hists_dict[names] = {'mean': None, 'std': None, 'bin_edges': None}
        FE_val[names] = {'mean': None, 'std': None}

    # These six dictionaries will be part of our output dict - They include overalls plus STD
    #    kde_gene_overall = {}
    #    kde_gene_std = {}
    #    kde_system_overall = {}
    #    kde_system_std = {}
    #    FE_curve_overall = {}
    #    FE_curve_std = {}
    #    hists_overall = {}
    #    hists_std = {}
    #    for names in enzymes_names:
    #        kde_gene_overall[names] = None
    #        kde_gene_std[names] = None
    #        kde_system_overall[names] = None
    ##       kde_system_std[names] = None
    #       FE_curve_overall[names] = None
    #       FE_curve_std[names] = None
    #       hists_overall[names] = None
    #       hists_std[names] = None

    # Let's prepare some data before running the parallelization
    # ------------------------------------------------------------------------------
    # Let's load the circuit, so we can extract the information that we need in an automatic way
    my_circuit = load_circuit(g_dict)

    # Get target site
    target_gene = [site for site in my_circuit.site_list if site.name == target_dict['target_gene']][0]
    RNAP_env = [environment for environment in my_circuit.environmental_list if environment.name == 'RNAP'][0]

    # Define x-axes
    x_system = get_interpolated_x(1, my_circuit.size)
    x_gene = get_interpolated_x(target_gene.start - RNAP_env.size, target_gene.end)

    # We need a list of items, so the pool can pass each item to the function
    Items = []
    for simulation_number in range(n_simulations):
        g_dict['n_simulations'] = simulation_number
        Item = {'global_conditions': g_dict, 'variations': variations_list[0]}  #because is list 0
        Items.append(Item)

    # TODO: AQUIMEQUDE
    #  1.- Set up loop.  - DONE
    #  2.- Write function sorts pool and gets: histograms, kdes, global, FE curve, FE, and CO - DONE
    #  3.- Write process that gets overalls and STDs - DONE
    #  4.- Test -DONE

    # Prepare dataframes to process parallel results
    # --------------------------------------------------------------
    # NOTE: We will collect all results in form of dataframes, then we will calculate the averages and STDs
    hist_df = {}
    kde_gene_df = {}
    kde_system_df = {}
    FE_curve_df = {}
    mean_FE_df = {}
    for p, name in enumerate(enzymes_names):
        hist_df[name] = None
        kde_gene_df[name] = None
        kde_system_df[name] = None
        FE_curve_df[name] = None
        mean_FE_df[name] = None

    # Run simulation and overalls
    # --------------------------------------------------------------
    # Create a multiprocessing pool
    pool = multiprocessing.Pool()

    # Run simulations for each n_set
    for n_set in range(n_sets):
        # Run simulations in parallel
        pool_results = pool.map(pt.single_simulation_w_variations_return_dfs, Items)

        # Process and extract results for the set
        nsupercoiling, nhist, nkde_gene, nkde_system, nFE_curve, nFE_val, ncorrelation = (
            extract_RNAPTrackingInfo_from_pool(
                pool_results, my_circuit, target_dict, reference_dict, x_gene, x_system))

        # Extract Global superhelical density
        # ------------------------------------------------------------------------------
        s_df = pd.DataFrame({'set_' + str(n_set): nsupercoiling})
        if n_set == 0:
            superhelical_df = s_df.copy()
        else:
            superhelical_df = pd.concat([superhelical_df, s_df], axis=1).reset_index(drop=True)

        # Extract correlation
        # ------------------------------------------------------------------------------
        c_df = pd.DataFrame({'set_' + str(n_set): ncorrelation})
        if n_set == 0:
            correlation_df = c_df.copy()
        else:
            correlation_df = pd.concat([correlation_df, c_df], axis=1).reset_index(drop=True)

        # Extract rest of data
        # ------------------------------------------------------------------------------
        for i, name in enumerate(enzymes_names):

            # Put data in data frame
            h_df = pd.DataFrame({'set_' + str(n_set): nhist[name]['count']})  # histogram counts
            k_g_df = pd.DataFrame({'set_' + str(n_set): nkde_gene[name]})  # kde_gene
            k_s_df = pd.DataFrame({'set_' + str(n_set): nkde_system[name]})   #kde_system
            FE_df = pd.DataFrame({'set_' + str(n_set): nFE_curve[name]})  # FE curve
            aFE_df = pd.DataFrame({'set_' + str(n_set): nFE_curve[name]})  # FE average value


            # Append data frame to a bigger data frame that we will use to process info
            if n_set == 0:
                hist_df[name] = h_df.copy()
                kde_gene_df[name] = k_g_df.copy()

                if name != 'RNAP':
                    kde_system_df[name] = k_s_df.copy()
                    FE_curve_df[name] = FE_df.copy()
                    mean_FE_df[name] = aFE_df.copy()
            else:
                hist_df = pd.concat([hist_df, h_df], axis=1).reset_index(drop=True)
                kde_gene_df = pd.concat([kde_gene_df, k_g_df], axis=1).reset_index(drop=True)

                if name != 'RNAP':
                    kde_system_df = pd.concat([kde_system_df, k_s_df], axis=1).reset_index(drop=True)
                    FE_curve_df = pd.concat([FE_curve_df, FE_df], axis=1).reset_index(drop=True)
                    mean_FE_df = pd.concat([mean_FE_df, aFE_df], axis=1).reset_index(drop=True)

    pool.close()  # Close the pool

    # Process data from the nsets and calculate overalls
    # --------------------------------------------------------------
    for i, name in enumerate(enzymes_names):

        # All enzymes/environmentals have these
        hists_dict[name]['mean'] = hist_df[name].mean(axis=1).to_numpy()
        hists_dict[name]['std'] = hist_df[name].std(axis=1).to_numpy()
        hists_dict[name]['bin_edges'] = h_df[name]['bin_edges']
        kde_gene[name]['mean'] = kde_gene_df[name].mean(axis=1).to_numpy()
        kde_gene[name]['std'] = kde_gene_df[name].std(axis=1).to_numpy()

        # RNAPs do not have fold-enrichment, only topos
        if name != 'RNAP':
            kde_system[name]['mean'] = kde_system_df[name].mean(axis=1).to_numpy()
            kde_system[name]['std'] = kde_system_df[name].std(axis=1).to_numpy()
            FE_curve[name]['mean'] = FE_curve_df[name].mean(axis=1).to_numpy()
            FE_curve[name]['std'] = FE_curve_df[name].std(axis=1).to_numpy()
            FE_val[name]['mean'] = mean_FE_df[name].mean(axis=1).to_numpy()
            FE_val[name]['std'] = mean_FE_df[name].std(axis=1).to_numpy()

    # Global supercoiling
    global_supercoiling['mean'] = superhelical_df.mean(axis=1).to_numpy()
    global_supercoiling['std'] = superhelical_df.std(axis=1).to_numpy()

    # For the mean correlation between sets
    correlation_sets['mean'] = correlation_df.mean(axis=1).to_numpy()
    correlation_sets['std'] = correlation_df.std(axis=1).to_numpy()

    # The overall correlation, calculated from the mean kdes
    buf_size = 10  # Used in calculating correlation. Is the number of data points ignored at the ends
    correlation_matrix = (
        np.corrcoef(kde_gene['topoI']['mean'][buf_size:-buf_size], kde_gene['RNAP']['mean'][buf_size:-buf_size]))
    overall_correlation = correlation_matrix[0,1]


    # Objective function
    #--------------------------------------------------------------------

    # Retrieve the target correlation and target FE
    target_CO = target_dict['target_CO']
    target_FE = target_dict['target_FE']

    my_objective = (target_FE - FE_val['topoI']['mean']) ** 2 + (target_CO - overall_correlation) ** 2

    # And prepare output dict
    #--------------------------------------------------------------------
    results_dict = {'kde_gene': kde_gene, 'kde_system': kde_system, 'FE_curve': FE_curve, 'hists': hists_dict,
                    'superhelical': global_supercoiling, 'FE_val': FE_val, 'correlation': correlation_sets}
    # These outputs dict contains the values used in the calculation of the objective function, plus all the
    # additional results.
    output_dict = {'FE': FE_val['topoI']['mean'], 'overall_correlation': overall_correlation, 'results': results_dict}

    return my_objective, output_dict


# Extracts information from parallelization that returns dfs, so it can be used in calibrating Topo I RNAPTracking
# pool_results = results from parallelization using pool.map
# Returns:
#     kde_gene = dict of KDEs interpolated to the target gene region
#     kde_system = dict of KDEs interpolated to the whole system
#     FE_dict = dict of fold-enrichment (FE) curves. RNAP does not have one
#     aFE_dict = dict of averaged fold-enrichment in the gene region.
#     hists_dict = dict of histograms
def extract_RNAPTrackingInfo_from_pool(pool_results, my_circuit, target_dict, reference_dict, x_gene, x_system):

    # Prepare output dicts
    # ------------------------------------------------------------------------------

    # These three dictionaries will be part of our output dict
    kde_gene = {}
    kde_system = {}
    FE_dict = {}  # curve
    aFE_dict = {}  # averaged fold-enrichment in the gene region.
    hists_dict = {}
    for names in target_dict['enzymes_names']:
        kde_gene[names] = None
        kde_system[names] = None
        FE_dict[names] = None
        aFE_dict[names] = None
        hists_dict[names] = None

    # Extract Global superhelical density
    # ------------------------------------------------------------------------------
    # Extract info from dataframes
    for j, out_dict in enumerate(pool_results):

        sites_df = out_dict['sites_df']
        mask = sites_df['type'] == 'circuit'
        superhelical = sites_df[mask]['superhelical'].to_numpy()

        simulation_df = pd.DataFrame({'simulation_' + str(j): superhelical})

        if j == 0:
            superhelical_output_df = simulation_df.copy()
        else:
            superhelical_output_df = pd.concat([superhelical_output_df, simulation_df], axis=1).reset_index(drop=True)

    superhelical_mean = superhelical_output_df.mean(axis=1).to_numpy()

    # Extract positions
    # ------------------------------------------------------------------------------

    # And filter results starting by name
    for i, name in enumerate(target_dict['enzymes_names']):

        all_positions = np.empty([1])

        # Extract info from dataframes
        for j, out_dict in enumerate(pool_results):
            enzymes_df = out_dict['enzymes_df']

            # Filter superhelical density
            mask = enzymes_df['name'] == name

            # Position
            # ------------------------------------------------------------------------------
            position = enzymes_df[mask]['position'].to_numpy()

            all_positions = np.concatenate([position, all_positions])

        # number of bins for histogram - which we will not plot
        nbins = calculate_number_nbins(my_circuit, name)

        # Calculate histogram
        counts, bin_edges = calculate_histogram(data=all_positions, nbins=nbins)
        hists_dict[name] = {'counts': counts, 'bin_edges': bin_edges}

        # print(name, len(all_positions))
        # Calculate KDE
        kde_x, kde_y = calculate_KDE(data=all_positions, nbins=nbins, scaled=True)

        # Let's interpolate for the gene - For every enzyme, let's interpolate signal in gene region for comparison
        kde_gene[name] = get_interpolated_kde(kde_x, kde_y, x_gene)

        # kde for the system (ony for topos)
        if name != 'RNAP':
            kde_system[name] = get_interpolated_kde(kde_x, kde_y, x_system)

            # Compute the fold-enrichment since we are here (no FE for RNAP)
            FE_dict[name] = kde_system[name] / reference_dict[0][name]

            # Compute average fold-enrichment FE
            aFE_dict[name] = calculate_mean_y_a_b(FE_dict[name], x_system, x_gene)

    # Compute correlation between topoI and RNAP signals along the gene
    # --------------------------------------------------------------------
    buf_size = 10  # Used in calculating correlation. Is the number of data points ignored at the ends
    correlation_matrix = np.corrcoef(kde_gene['topoI'][buf_size:-buf_size], kde_gene['RNAP'][buf_size:-buf_size])
    correlation = correlation_matrix[0, 1]

    return superhelical_mean, hists_dict, kde_gene, kde_system, FE_dict, aFE_dict, correlation

# This calculates the mean of y within the ranges z[0] and z[-1] but y is paired with x in the form of [x,y].
# We use this function to calculate the Fold-enrichment of topo I within the gene region specified by the ends of z.
def calculate_mean_y_a_b(y, x, z):

    # Find the index of the value in x that is closest to z[0]
    a = np.abs(x - z[0]).argmin()

    # Find the index of the value in x that is closest to z[-1]
    b = np.abs(x - z[-1]).argmin()

    # Ensure a is less than or equal to b
    if a > b:
        a, b = b, a

    # Calculate the mean of y in the range [a, b]
    mean_y = np.mean(y[a:b + 1])
    return mean_y

# This function just loads the circuit from the global dict and returns the circuit (without running it).
# I put as a function so in case you just need to process the circuit to extract some global data, you don't
# actually need to write the whole thing. So it is just to save some space and clean up a bit the functions.
def load_circuit(global_dict):
    my_circuit = Circuit(circuit_filename=global_dict['circuit_filename'], sites_filename=global_dict['sites_filename'],
                         enzymes_filename=global_dict['enzymes_filename'],
                         environment_filename=global_dict['environment_filename'],
                         output_prefix=global_dict['output_prefix'], frames=global_dict['frames'],
                         series=global_dict['series'], continuation=global_dict['continuation'],
                         dt=global_dict['dt'])

    return my_circuit


# Calculates the number of bins for calculating histogram. I put as a function so we can use the same criteria
# for all cases.
# my_circuit is the circuit, and name is the environmental name, e.g., topoI, gyrase, RNAP
def calculate_number_nbins(my_circuit, name):
    factor = .5
    gene_factor = .5  #.4#.2#.35
    # Get environmental object
    my_environmental = [environmental for environmental in my_circuit.environmental_list if environmental.name == name][
        0]

    # So, if topos or environmentals that bind the whole DNA
    if 'DNA_' in my_environmental.site_type:
        nbins = int(
            factor * my_circuit.size / my_environmental.size)  # number of bins. - note that this only applies for topos
    else:  # For genes, etc...
        #size = abs(my_environmental.site_list[0].start - my_environmental.site_list[0].end)
        size = abs(my_environmental.site_list[0].start - my_environmental.size - my_environmental.site_list[0].end)
        nbins = int(gene_factor * size / my_environmental.size)  # because genes are smaller

    return nbins


def calculate_histogram(data, nbins):
    # Calculate the histogram without normalization
    counts, bin_edges = np.histogram(data, bins=nbins, density=False)

    return counts, bin_edges


# Calculates kde from a histogram constructed for data with nbins
def calculate_KDE(data, nbins, scaled=True):
    # Calculate the histogram without normalization
    counts, bin_edges = np.histogram(data, bins=nbins, density=False)

    # Calculate the bin width
    bin_width = bin_edges[1] - bin_edges[0]

    # calculate kde with scipy
    kde = gaussian_kde(data)
    kde_x = np.linspace(min(data), max(data), nbins)
    kde_y = kde(kde_x)

    # Scale KDE values
    if scaled:
        kde_y_scaled = kde_y * len(data) * bin_width
        kde_y = kde_y_scaled

    return kde_x, kde_y


def get_interpolated_kde(kde_x, kde_y, x_interpolated):
    # Create interpolation function
    interp_fun = interp1d(kde_x, kde_y, kind='linear', fill_value='extrapolate')  # everything default

    # Get interpolated y-values
    y_interpolated = interp_fun(x_interpolated)
    return y_interpolated


# This one does not actually interpolates the x-axis, but it get's the axis for which the kde_y will be interpolated to
def get_interpolated_x(start, end, x_spacing=10):
    # Define x-axes
    x_interpolated = np.arange(start, end, x_spacing)

    return x_interpolated
