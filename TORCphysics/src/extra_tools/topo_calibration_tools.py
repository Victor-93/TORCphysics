import numpy as np
import multiprocessing
import multiprocessing.pool
from TORCphysics import parallelization_tools as pt
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
from TORCphysics import Circuit
from TORCphysics import analysis as ann
import pandas as pd


# Stuff added to make parallel / parallel processes
# -----------------------------------------------------------------------------------
class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, val):
        pass


#class NoDaemonProcessPool(multiprocessing.pool.Pool):
class MyPool(multiprocessing.pool.Pool):

    def Process(self, *args, **kwds):
        #proc = super(NoDaemonProcessPool, self).Process(*args, **kwds)
        proc = super(MyPool, self).Process(*args, **kwds)
        proc.__class__ = NoDaemonProcess

        return proc


#class NoDaemonProcess(multiprocessing.Process):
#    # make 'daemon' attribute always return False
#    def _get_daemon(self):
#        return False
#    def _set_daemon(self, value):
#        pass
#    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
#class MyPool(multiprocessing.pool.Pool):
#    Process = NoDaemonProcess


# Custom Pool class to allow nested parallelism
#class NoDaemonProcess(multiprocessing.Process):
#    # Make 'daemon' attribute always return False
#    @property
#    def daemon(self):
#        return False#
#
#    @daemon.setter
#    def daemon(self, value):
#        pass

#class NoDaemonPool(multiprocessing.pool.Pool):
#    Process = NoDaemonProcess

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
        c_df = pd.DataFrame({'set_' + str(n_set): [ncorrelation]})
        if n_set == 0:
            correlation_df = c_df.copy()
        else:
            correlation_df = pd.concat([correlation_df, c_df], axis=1).reset_index(drop=True)

        # Extract rest of data
        # ------------------------------------------------------------------------------
        for i, name in enumerate(enzymes_names):

            # Put data in data frame
            h_df = pd.DataFrame({'set_' + str(n_set): nhist[name]['counts']})  # histogram counts
            k_g_df = pd.DataFrame({'set_' + str(n_set): nkde_gene[name]})  # kde_gene
            if name != 'RNAP':
                k_s_df = pd.DataFrame({'set_' + str(n_set): nkde_system[name]})  #kde_system
                FE_df = pd.DataFrame({'set_' + str(n_set): nFE_curve[name]})  # FE curve
                aFE_df = pd.DataFrame({'set_' + str(n_set): [nFE_val[name]]})  # FE average value

            # Append data frame to a bigger data frame that we will use to process info
            if n_set == 0:
                hist_df[name] = h_df.copy()
                kde_gene_df[name] = k_g_df.copy()

                hists_dict[name]['bin_edges'] = nhist[name]['bin_edges']  #Since we are here, let's do it now.

                if name != 'RNAP':
                    kde_system_df[name] = k_s_df.copy()
                    FE_curve_df[name] = FE_df.copy()
                    mean_FE_df[name] = aFE_df.copy()
            else:
                hist_df[name] = pd.concat([hist_df[name], h_df], axis=1).reset_index(drop=True)
                kde_gene_df[name] = pd.concat([kde_gene_df[name], k_g_df], axis=1).reset_index(drop=True)

                if name != 'RNAP':
                    kde_system_df[name] = pd.concat([kde_system_df[name], k_s_df], axis=1).reset_index(drop=True)
                    FE_curve_df[name] = pd.concat([FE_curve_df[name], FE_df], axis=1).reset_index(drop=True)
                    mean_FE_df[name] = pd.concat([mean_FE_df[name], aFE_df], axis=1).reset_index(drop=True)

    pool.close()  # Close the pool

    # Process data from the nsets and calculate overalls
    # --------------------------------------------------------------
    for i, name in enumerate(enzymes_names):

        # All enzymes/environmentals have these
        hists_dict[name]['mean'] = hist_df[name].mean(axis=1).to_numpy()
        hists_dict[name]['std'] = hist_df[name].std(axis=1).to_numpy()
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
    overall_correlation = correlation_matrix[0, 1]

    # Let's calculate the correlation between the RNAP KDE and the reference RNA (experimental fitted to the TU)
    correlation_matrix = np.corrcoef(kde_gene['RNAP']['mean'], reference_dict[0]['RNAP'])
    RNAP_correlation = correlation_matrix[0, 1]

    # Objective function
    #--------------------------------------------------------------------

    # Retrieve the target correlation and target FE
    target_CO = target_dict['target_CO']
    target_FE = target_dict['target_FE']

    my_objective = (target_FE - FE_val['topoI']['mean']) ** 2 + (target_CO - overall_correlation) ** 2
    my_objective = my_objective[0]

    # And prepare output dict
    #--------------------------------------------------------------------
    results_dict = {'kde_gene': kde_gene, 'kde_system': kde_system, 'FE_curve': FE_curve, 'hists_dict': hists_dict,
                    'superhelical_dict': global_supercoiling, 'FE_val': FE_val, 'correlation': correlation_sets,
                    'RNAP_correlation': RNAP_correlation}

    # These outputs dict contains the values used in the calculation of the objective function, plus all the
    # additional results.
    output_dict = {'FE': FE_val['topoI']['mean'], 'overall_correlation': overall_correlation, 'results': results_dict,
                   'objective': my_objective}  # Let's also include the objective so we have all the info available.

    return my_objective, output_dict


# This function is very similar than single_case_RNAPTracking_calibration_nsets_2scheme().
# It also performs multiple sets of simulations to smooth the KDEs, but this time it uses a different
# parallelization scheme (that's why is called *_2scheme).
# Essentially, the scheme consists in creating two different parallelization Pools(), one outer and one inner.
# In the outter pool, the number of n_sets are partitioned, and the number of n_workers are distributed to these
# outer pools. Within each pool, an inner pool is created where it runs n_inner_workers simulations. This is done
# n_subset times. So, in total, we run n_sets*n_subsets number of sets to smooth the KDEs.
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
# WARNING: It is very important that the number of n_workers is within the capabilities of your system.
def single_case_RNAPTracking_calibration_nsets_2scheme(global_dict_list, variations_list, reference_dict, target_dict):
    g_dict = dict(global_dict_list[0])  # Just one system
    n_workers = g_dict['n_workers']  # Total number of workers (cpus)
    n_sets = g_dict['n_sets']  # Number of outer sets
    n_subsets = g_dict['n_subsets']  # Number of simulations per set
    n_inner_workers = g_dict['n_inner_workers']  # Number of workers per inner pool

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

    # Prepare dataframes to process parallel results
    # --------------------------------------------------------------
    # NOTE: We will collect all results in form of dataframes, then we will calculate the averages and STDs
    hist_df = {}
    kde_gene_df = {}
    kde_system_df = {}
    FE_curve_df = {}
    mean_FE_df = {}
    bin_edges = {}
    for p, name in enumerate(enzymes_names):
        hist_df[name] = None
        kde_gene_df[name] = None
        kde_system_df[name] = None
        FE_curve_df[name] = None
        mean_FE_df[name] = None
        bin_edges[name] = None

    # Prepare items for parallelization
    # --------------------------------------------------------------
    # We need a list of items, so the pool can pass each item to the function
    s = 0
    Item_set = []
    for n_set in range(n_sets):
        Item_subset = []
        for n_subset in range(n_subsets):
            Items = []
            for n_inner_worker in range(n_inner_workers):
                g_dict['n_simulations'] = s  # This is the simulation number
                Item = {'global_conditions': g_dict.copy(), 'variations': variations_list[0]}
                Items.append(Item)
                s += 1
            Item_subset.append(Items)
        Item_set.append(Item_subset)

    processing_info_dict = {'circuit': my_circuit, 'target_dict': target_dict, 'reference_dict': reference_dict,
                            'x_gene': x_gene, 'x_system': x_system}
    # Launch parellelization scheme v2
    # --------------------------------------------------------------

    # Create a multiprocessing pool for the outer loop
    pool = MyPool(n_sets)
    Item_forpool = [(n_set, n_subsets, Item_set[n_set], n_inner_workers, processing_info_dict) for n_set in
                    range(n_sets)]

    results = pool.map(
        process_set, Item_forpool)

    pool.close()
    #with NoDaemonPool(n_sets) as outer_pool:
    # with multiprocessing.Pool(n_sets) as outer_pool:
    # Run the outer loop in parallel
    #    results = outer_pool.starmap(
    #        process_set, [(n_set, n_subsets, Item_set, n_inner_workers, processing_info_dict) for n_set in range(n_sets)])

    # Let's retrieve the results and organise them
    # --------------------------------------------------------------
    for j, processed_dict in enumerate(results):

        if j == 0:
            for i, name in enumerate(enzymes_names):
                hist_df[name] = processed_dict['hist_df'][name].copy()
                kde_gene_df[name] = processed_dict['kde_gene_df'][name].copy()
                hists_dict[name]['bin_edges'] = processed_dict['bin_edges'][name]
                if name != 'RNAP':
                    kde_system_df[name] = processed_dict['kde_system_df'][name].copy()
                    FE_curve_df[name] = processed_dict['FE_curve_df'][name].copy()
                    mean_FE_df[name] = processed_dict['mean_FE_df'][name].copy()
            superhelical_df = processed_dict['superhelical_df'].copy()
            correlation_df = processed_dict['correlation_df'].copy()
        else:
            for i, name in enumerate(enzymes_names):
                hist_df[name] = pd.concat([hist_df[name], processed_dict['hist_df'][name]], axis=1).reset_index(
                    drop=True)
                kde_gene_df[name] = pd.concat([kde_gene_df[name], processed_dict['kde_gene_df'][name]],
                                              axis=1).reset_index(drop=True)
                if name != 'RNAP':
                    kde_system_df[name] = pd.concat([kde_system_df[name], processed_dict['kde_system_df'][name]],
                                                    axis=1).reset_index(drop=True)
                    FE_curve_df[name] = pd.concat([FE_curve_df[name], processed_dict['FE_curve_df'][name]],
                                                  axis=1).reset_index(drop=True)
                    mean_FE_df[name] = pd.concat([mean_FE_df[name], processed_dict['mean_FE_df'][name]],
                                                 axis=1).reset_index(drop=True)
            superhelical_df = pd.concat([superhelical_df, processed_dict['superhelical_df']], axis=1).reset_index(
                drop=True)
            correlation_df = pd.concat([correlation_df, processed_dict['correlation_df']], axis=1).reset_index(
                drop=True)

    # Process data from the nsets and calculate overalls
    # --------------------------------------------------------------
    for i, name in enumerate(enzymes_names):

        # All enzymes/environmentals have these
        hists_dict[name]['mean'] = hist_df[name].mean(axis=1).to_numpy()
        hists_dict[name]['std'] = hist_df[name].std(axis=1).to_numpy()
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
    overall_correlation = correlation_matrix[0, 1]

    # Let's calculate the correlation between the RNAP KDE and the reference RNA (experimental fitted to the TU)
    correlation_matrix = np.corrcoef(kde_gene['RNAP']['mean'], reference_dict[0]['RNAP'])
    RNAP_correlation = correlation_matrix[0, 1]

    # Objective function
    #--------------------------------------------------------------------

    # Retrieve the target correlation and target FE
    target_CO = target_dict['target_CO']
    target_FE = target_dict['target_FE']

    my_objective = (target_FE - FE_val['topoI']['mean']) ** 2 + (target_CO - overall_correlation) ** 2
    my_objective = my_objective[0]

    # And prepare output dict
    #--------------------------------------------------------------------
    results_dict = {'kde_gene': kde_gene, 'kde_system': kde_system, 'FE_curve': FE_curve, 'hists_dict': hists_dict,
                    'superhelical_dict': global_supercoiling, 'FE_val': FE_val, 'correlation': correlation_sets,
                    'RNAP_correlation': RNAP_correlation}

    # These outputs dict contains the values used in the calculation of the objective function, plus all the
    # additional results.
    output_dict = {'FE': FE_val['topoI']['mean'], 'overall_correlation': overall_correlation, 'results': results_dict,
                   'objective': my_objective}  # Let's also include the objective so we have all the info available.

    return my_objective, output_dict


def process_set(item_pool):
    #    n_set, n_subsets, Items_subset, n_inner_workers, processing_info_dict
    n_set = item_pool[0]
    n_subsets = item_pool[1]
    Items_subset = item_pool[2]
    n_inner_workers = item_pool[3]
    processing_info_dict = item_pool[4]

    # Extract info we need for processing
    # --------------------------------------------------------------
    my_circuit = processing_info_dict['circuit']
    target_dict = processing_info_dict['target_dict']
    reference_dict = processing_info_dict['reference_dict']
    x_gene = processing_info_dict['x_gene']
    x_system = processing_info_dict['x_system']
    enzymes_names = target_dict['enzymes_names']

    # Prepare dataframes to process parallel results
    # --------------------------------------------------------------
    # NOTE: We will collect all results in form of dataframes, then we will calculate the averages and STDs
    hist_df = {}
    kde_gene_df = {}
    kde_system_df = {}
    FE_curve_df = {}
    mean_FE_df = {}
    bin_edges = {}
    for p, name in enumerate(enzymes_names):
        hist_df[name] = None
        kde_gene_df[name] = None
        kde_system_df[name] = None
        FE_curve_df[name] = None
        mean_FE_df[name] = None
        bin_edges[name] = None

    # Create a multiprocessing pool for the inner loop
    # with multiprocessing.Pool(n_inner_workers) as inner_pool:
    #    with NoDaemonPool(n_inner_workers) as inner_pool:

    inner_pool = multiprocessing.Pool(n_inner_workers)

    for n in range(n_subsets):

        Items = Items_subset[n]
        # Run simulations in parallel within this subset
        pool_results = inner_pool.map(pt.single_simulation_w_variations_return_dfs, Items)

        # Process and extract results for the set
        nsupercoiling, nhist, nkde_gene, nkde_system, nFE_curve, nFE_val, ncorrelation = (
            extract_RNAPTrackingInfo_from_pool(
                pool_results, my_circuit, target_dict, reference_dict, x_gene, x_system))

        # Extract Global superhelical density
        # ------------------------------------------------------------------------------
        s_df = pd.DataFrame({'set_' + str(n): nsupercoiling})
        if n == 0:
            superhelical_df = s_df.copy()
        else:
            superhelical_df = pd.concat([superhelical_df, s_df], axis=1).reset_index(drop=True)

        # Extract correlation
        # ------------------------------------------------------------------------------
        c_df = pd.DataFrame({'set_' + str(n): [ncorrelation]})
        if n == 0:
            correlation_df = c_df.copy()
        else:
            correlation_df = pd.concat([correlation_df, c_df], axis=1).reset_index(drop=True)

        # Extract rest of data
        # ------------------------------------------------------------------------------
        for i, name in enumerate(enzymes_names):

            # Put data in data frame
            h_df = pd.DataFrame({'set_' + str(n): nhist[name]['counts']})  # histogram counts
            k_g_df = pd.DataFrame({'set_' + str(n): nkde_gene[name]})  # kde_gene
            if name != 'RNAP':
                k_s_df = pd.DataFrame({'set_' + str(n): nkde_system[name]})  # kde_system
                FE_df = pd.DataFrame({'set_' + str(n): nFE_curve[name]})  # FE curve
                aFE_df = pd.DataFrame({'set_' + str(n): [nFE_val[name]]})  # FE average value

            # Append data frame to a bigger data frame that we will use to process info
            if n == 0:
                hist_df[name] = h_df.copy()
                kde_gene_df[name] = k_g_df.copy()

                bin_edges[name] = nhist[name]['bin_edges']  # Since we are here, let's do it now.

                if name != 'RNAP':
                    kde_system_df[name] = k_s_df.copy()
                    FE_curve_df[name] = FE_df.copy()
                    mean_FE_df[name] = aFE_df.copy()
            else:
                hist_df[name] = pd.concat([hist_df[name], h_df], axis=1).reset_index(drop=True)
                kde_gene_df[name] = pd.concat([kde_gene_df[name], k_g_df], axis=1).reset_index(drop=True)

                if name != 'RNAP':
                    kde_system_df[name] = pd.concat([kde_system_df[name], k_s_df], axis=1).reset_index(drop=True)
                    FE_curve_df[name] = pd.concat([FE_curve_df[name], FE_df], axis=1).reset_index(drop=True)
                    mean_FE_df[name] = pd.concat([mean_FE_df[name], aFE_df], axis=1).reset_index(drop=True)

    processed = {'hist_df': hist_df, 'kde_gene_df': kde_gene_df, 'kde_system_df': kde_system_df,
                 'FE_curve_df': FE_curve_df, 'mean_FE_df': mean_FE_df, 'superhelical_df': superhelical_df,
                 'correlation_df': correlation_df, 'bin_edges': bin_edges}
    return processed


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


# This function processes the results from pool (pool_outputs), and given a list of site_names, it returns
#  a dict of dictionaries with each entry being a DF with average values of
#  superhelical, binding, unbinding and #enzymes.
#  Given a list of enzyme_names, it returns a dictionary with numpy arrays of their KDEs.
#  All these parameters are averaged values
#  Outputs: sites_dict{'superhelical', 'binding', 'unbinding', '#enzymes' }.
#           For each of the site_dataframes, columns have the form: site_name_mean, site_name_std; for the sites in
#           site_names list
#           enzymes_dict{'enzyme_names'} - each one containin a KDE  - Or maybe calculated after?
def process_pools_get_avgs(pool_outputs, site_names, enzyme_names, my_circuit):
    # Create a range of values for the 'frames' column
    frame_values = list(range(my_circuit.frames + 1))

    # Calculate the corresponding time values using dt
    time_values = [my_circuit.dt * frame for frame in frame_values]

    # Create a DataFrame with 'frames' and 'time' columns
    time_df = pd.DataFrame({'frames': frame_values, 'time': time_values})

    # Let's do it for sites first!
    # ------------------------------------------------------------------------------
    # Prepare output dataframe
    superhelical_output_df = time_df.copy()
    binding_output_df = time_df.copy()
    unbinding_output_df = time_df.copy()
    nenzymes_output_df = time_df.copy()

    # And filter site results starting by name
    for name in site_names:
        for n, out_dict in enumerate(pool_outputs):

            #Get output dfs
            sites_df = out_dict['sites_df']

            # Filter superhelical density
            if name == 'circuit':
                mask = sites_df['type'] == name
            else:
                mask = sites_df['name'] == name

            # Collect measurements
            # ------------------------------------------------------------------------------
            superhelical = sites_df[mask]['superhelical'].to_numpy()
            binding_event = sites_df[mask]['binding'].to_numpy()
            unbinding_event = sites_df[mask]['unbinding'].to_numpy()
            number_enzymes = sites_df[mask]['#enzymes'].to_numpy()

            # Make them DFs and sort by simulation
            superhelical = pd.DataFrame({'simulation_' + str(n): superhelical})
            binding_event = pd.DataFrame({'simulation_' + str(n): binding_event})
            unbinding_event = pd.DataFrame({'simulation_' + str(n): unbinding_event})
            number_enzymes = pd.DataFrame({'simulation_' + str(n): number_enzymes})

            if n == 0:
                superhelical_df = superhelical.copy()
                binding_event_df = binding_event.copy()
                unbinding_event_df = unbinding_event.copy()
                number_enzymes_df = number_enzymes.copy()
            else:
                superhelical_df = pd.concat([superhelical_df, superhelical], axis=1).reset_index(drop=True)
                binding_event_df = pd.concat([binding_event_df, binding_event], axis=1).reset_index(drop=True)
                unbinding_event_df = pd.concat([unbinding_event_df, unbinding_event], axis=1).reset_index(drop=True)
                number_enzymes_df = pd.concat([number_enzymes_df, number_enzymes], axis=1).reset_index(drop=True)

        # Calculate averages
        # ------------------------------------------------------------------------------
        #Superhelical
        mean_df = superhelical_df.mean(axis=1).to_frame(name=name + '_mean')
        std_df = superhelical_df.std(axis=1).to_frame(name=name + '_std')
        superhelical_output_df = pd.concat([superhelical_output_df, mean_df, std_df], axis=1)  #.reset_index(True)

        # Binding
        mean_df = binding_event_df.mean(axis=1).to_frame(name=name + '_mean')
        std_df = binding_event_df.std(axis=1).to_frame(name=name + '_std')
        binding_output_df = pd.concat([binding_output_df, mean_df, std_df], axis=1)  #.reset_index(True)

        # Unbinding
        mean_df = unbinding_event_df.mean(axis=1).to_frame(name=name + '_mean')
        std_df = unbinding_event_df.std(axis=1).to_frame(name=name + '_std')
        unbinding_output_df = pd.concat([unbinding_output_df, mean_df, std_df], axis=1)

        # #Enzymes
        mean_df = number_enzymes_df.mean(axis=1).to_frame(name=name + '_mean')
        std_df = number_enzymes_df.std(axis=1).to_frame(name=name + '_std')
        nenzymes_output_df = pd.concat([nenzymes_output_df, mean_df, std_df], axis=1)

    # And the sites_dict output dictionary!
    sites_dict = {'superhelical': superhelical_output_df, 'binding': binding_output_df,
                  'unbinding': unbinding_output_df, '#enzymes': nenzymes_output_df}

    # For the enzymes.
    # ------------------------------------------------------------------------------
    enzymes_dict = {}
    for name in enzyme_names:

        all_positions = np.empty([1])

        # Extract info from dataframes
        for n, out_dict in enumerate(pool_outputs):
            enzymes_df = out_dict['enzymes_df']

            # Filter
            mask = enzymes_df['name'] == name

            # Position
            # ------------------------------------------------------------------------------
            position = enzymes_df[mask]['position'].to_numpy()

            all_positions = np.concatenate([position, all_positions])

        # number of bins for histogram - which we will not plot
        nbins = calculate_number_nbins(my_circuit, name)

        # Calculate histogram
        counts, bin_edges = calculate_histogram(data=all_positions, nbins=nbins)

        # Calculate KDE
        kde_x, kde_y = calculate_KDE(data=all_positions, nbins=nbins, scaled=True)

        enzymes_dict[name + '_counts'] = counts
        enzymes_dict[name + '_bin_edges'] = bin_edges
        enzymes_dict[name + '_kde_x'] = kde_x
        enzymes_dict[name + '_kde_y'] = kde_y

    return sites_dict, enzymes_dict


# TODO: This function only calculates the susceptibilities given the results from Pool()
# The production rate is calculated as: prod_rate = sum(unbinding_events)/total_time
#  In other words, counts how many complete cycles or outputs were done on average per circuit in the total
#  simulation time.
# Inputs: sites_names = list with site names to calculate susceptibility from
# Returns: Dictionary with susceptibilities, in the form {site_name: susceptibility}
def process_pools_calculate_production_rate(pool_outputs, site_names, my_circuit):
    total_time = my_circuit.dt * my_circuit.frames

    time = np.arange(0, total_time + my_circuit.dt, my_circuit.dt)

    production_dict = {}  # This is the output

    # Calculation time!
    # ------------------------------------------------------------------------------

    # And filter site results starting by name
    for name in site_names:
        for n, out_dict in enumerate(pool_outputs):

            #Get output dfs
            sites_df = out_dict['sites_df']

            # Filter superhelical density
            if name == 'circuit':
                mask = sites_df['type'] == name
            else:
                mask = sites_df['name'] == name

            # Collect measurements
            # ------------------------------------------------------------------------------
            unbinding_event = sites_df[mask]['unbinding'].to_numpy()

            # Make them DFs and sort by simulation
            unbinding_event = pd.DataFrame({'simulation_' + str(n): unbinding_event})

            if n == 0:
                unbinding_event_df = unbinding_event.copy()
            else:
                unbinding_event_df = pd.concat([unbinding_event_df, unbinding_event], axis=1).reset_index(drop=True)

        # Calculate averages
        # ------------------------------------------------------------------------------
        mean = unbinding_event_df.mean(axis=1).to_numpy()

        sum_time, log_sum, prod_rate = ann.calculate_steady_state_initiation_curve(mean, time,
                                                                                   ta=int(my_circuit.frames / 2), tb=-1)

        #prod_rate = sum(mean)/total_time

        production_dict[name] = prod_rate

    return production_dict


# TODO: Do a bit of cleaning
def gene_architecture_calibration_nsets(big_global_list, big_variation_list, experimental_curves, parallel_info,
                                        additional_results=False):
    # Let's unpack the info
    ncases = len(big_global_list)
    n_sets = parallel_info['n_sets']
    n_subsets = parallel_info['n_subsets']
    n_inner_workers = parallel_info['n_inner_workers']
    objective = 0  # Total objective function
    output_list = []  # This list will be made out of a list of cases, each case is made of  a dict with
    # objective and data keys. In data, it gives a list of dicts, each one containing the distance, and the results

    # Go through list of cases
    # --------------------------------------------------------------
    # Create a multiprocessing pool for the outer loop
    pool = MyPool(n_sets)  # I do it now so I don't have to open it for every case

    # Go through each case
    for icase in range(ncases):

        global_list = big_global_list[icase]
        variation_list = big_variation_list[icase]
        exp_curve = experimental_curves[icase]

        distances = exp_curve['distance']

        # Prepare items for parallelization
        # --------------------------------------------------------------
        # We need a list of items, so the pool can pass each item to the function
        # Even though all inner workers within the same n_set have the same global and variation conditions,
        # We need each simulation to have their own ID (so they also have their own random seed)
        s = 0
        Item_set = []
        processing_dict_list = []  # This dict will help us process the outputs
        for j in range(len(distances)):
            g_dict = global_list[j]
            Item_subset = []
            for n_subset in range(n_subsets):
                Items = []
                for n_inner_worker in range(n_inner_workers):
                    g_dict['n_simulations'] = s  # This is the simulation number
                    Item = {'global_conditions': g_dict.copy(), 'variations': variation_list.copy()}
                    Items.append(Item)
                    s += 1
                Item_subset.append(Items)
            Item_set.append(Item_subset)

            # Sort the processing info
            # --------------------------------------------------
            my_circuit = load_circuit(g_dict)

            # Get target site
            target_gene = [site for site in my_circuit.site_list if site.name == 'reporter'][0]
            RNAP_env = [environment for environment in my_circuit.environmental_list if environment.name == 'RNAP'][0]

            # Define x-axes
            x_system = get_interpolated_x(1, my_circuit.size)
            x_gene = get_interpolated_x(target_gene.start - RNAP_env.size, target_gene.end)

            p_dict = {'circuit': my_circuit, 'x_gene': x_gene, 'x_system': x_system, 'site_name': 'reporter',
                      'additional_results': additional_results}
            processing_dict_list.append(p_dict)

        # Launch parellelization scheme
        # --------------------------------------------------------------
        Item_forpool = [(n_subsets, n_inner_workers, Item_set[j], processing_dict_list[j]) for j in
                        range(len(distances))]

        results = pool.map(
            gene_architecture_run_pool, Item_forpool)

        # Process and calculate objective
        # --------------------------------------------------------------
        prod_rates = [result['prod_rate'][0] for result in results]
        susceptibility = prod_rates / prod_rates[4]

        my_objective = np.sum((exp_curve['Signal'] - susceptibility) ** 2)  # This is the objective of the icase
        objective += my_objective  # This is the total objective

        # Prepare output - We will just add the distance
        out_list = []  # output for the icase
        for i, result in enumerate(results):
            out_dict = {}
            out_dict['distance'] = distances[i]
            out_dict['result'] = result
            #out_dict['objective'] = my_objective  # Save the objective of the current case
            out_list.append(out_dict)
        # output_list.append(out_list)
        output_list.append({'objective': my_objective, 'data': out_list})

    # Close pool
    # --------------------------------------------------------------
    pool.close()

    return objective, output_list


def gene_architecture_run_pool(item_pool):
    n_subsets = item_pool[0]  # This is how many sets of simulations are run
    n_inner_workers = item_pool[1]  # Number of workers
    Items_subset = item_pool[2]  # Conditions
    processing_dict = item_pool[3]
    additional_results = processing_dict['additional_results']

    # Our output will be a dictionary
    output = {}

    # Total number of simulations ran are n_subsets * n_inner_workers

    # Run inner parallelization
    # ---------------------------------------------------------
    # Create a multiprocessing pool for the inner loop
    inner_pool = multiprocessing.Pool(n_inner_workers)

    pool_list = []  # This collects all the results.
    for n in range(n_subsets):
        Items = Items_subset[n]
        # Run simulations in parallel within this subset
        pool_results = inner_pool.map(pt.single_simulation_w_variations_return_dfs, Items)

        pool_list += pool_results

    processing_items = []
    for pool_results in pool_list:
        processing_items.append((pool_results, processing_dict))

    # Process results with another parallelization process
    # ---------------------------------------------------------
    # Let's calculate one of each, and do the average and get std
    output_pools = inner_pool.map(gene_architecture_process_pool, processing_items)

    # Get averages and stds
    # ---------------------------------------------------------

    # Production rate -------------------------
    # Extract prod_rate values from the output_pools list
    prod_rates = [result[0] for result in output_pools]

    # Convert to a NumPy array
    prod_rates_array = np.array(prod_rates)

    # Calculate mean and standard deviation
    mean_prod_rate = np.mean(prod_rates_array)
    std_prod_rate = np.std(prod_rates_array)

    prod_rate = np.array([mean_prod_rate, std_prod_rate])

    # Additional results -------------------------
    if additional_results:
        ares = [result[1] for result in output_pools]  # The dict returned by gene_architecture_process_pool

        # Procesing DFs and calculating mean,std
        for key in ['global_superhelical', 'local_superhelical', 'binding', 'unbinding']:
            # Extract the NumPy arrays from each dictionary
            arrays = [d[key] for d in ares]

            # Stack the arrays into a 2D NumPy array
            stacked_arrays = np.vstack(arrays)

            # Calculate the mean and standard deviation along the first axis
            mean_array = np.mean(stacked_arrays, axis=0)
            std_array = np.std(stacked_arrays, axis=0)

            # And add it to the dict
            output[key] = np.array([mean_array, std_array])

        # Now for the KDEs
        kde_dict = {}
        arrays = [d['KDE'] for d in ares]  # We have all the KDEs
        for name in ['RNAP', 'topoI', 'gyrase']:
            # Extract KDE (kde_y) from the molecule with "name"
            arrays_kde_y = [d[name]['kde_y'] for d in arrays]  # We collects all KDE_y of name

            # Stack the arrays into a 2D NumPy array
            stacked_arrays = np.vstack(arrays_kde_y)

            # Calculate the mean and standard deviation along the first axis
            mean_array = np.mean(stacked_arrays, axis=0)
            std_array = np.std(stacked_arrays, axis=0)

            kde_y = np.array([mean_array, std_array])
            kde_x = arrays[0][name]['kde_x']

            # And add it to the dict
            kde_dict[name] = {'kde_y': kde_y, 'kde_x': kde_x}

        # And add it to our
        output['KDE'] = kde_dict

    output['prod_rate'] = prod_rate
    # process_pools_genearch(pool_list, processing_dict)

    # In the meantime, let's just return the production rate
    return output


# This one processes the outputs of the pool (used in the gene architecture experiment).
# Returns: production rate (prod_rate), and accordition additional_results, it get's the following outputs:
# If additional_results: False, we only return production rates (prod_rate), defined as the number of
# transcripts divided by the elapsed time
# If additional_results: True, returns KDEs of RNAP, topoI and gyrase. It also returns global superhelical level,
# and local superhelical level at site, as well as the binding and unbinding arrays at the site
def gene_architecture_process_pool(item_pool):
    # Unpack item_pool
    pool_result = item_pool[0]  # These are the df's produced by TORCphysics
    processing_dict = item_pool[1]

    # Unpack processing dict
    additional_results = processing_dict['additional_results']
    my_circuit = processing_dict['circuit']
    x_gene = processing_dict['x_gene']
    x_system = processing_dict['x_system']
    site_name = processing_dict['site_name']

    # Let's calculate the production rate
    # -----------------------------------------------------------------------------------------------------------------
    total_time = my_circuit.dt * my_circuit.frames

    # Ranges used for calculating unbinding events
    a = 0  #int(my_circuit.frames / 2)
    b = -1

    # Get unbinding event
    sites_df = pool_result['sites_df']
    mask = sites_df['name'] == site_name
    unbinding_event = sites_df[mask]['unbinding'].to_numpy()

    # Calculate production rate (prod_rate)
    prod_rate = np.sum(unbinding_event[a:b]) / total_time  #/ 2
    #production_rate_dict = process_pools_calculate_production_rate([pool_result], [site_name], my_circuit)
    #production_rate = production_rate_dict[site_name]

    additional_dict = None
    # Let's do the processing
    if additional_results:

        additional_dict = {}
        # Collect measurements
        # ------------------------------------------------------------------------------
        # Site measurements
        local_superhelical = sites_df[mask]['superhelical'].to_numpy()
        binding_event = sites_df[mask]['binding'].to_numpy()

        # Global measurements
        mask_circuit = sites_df['type'] == 'circuit'
        global_superhelical = sites_df[mask_circuit]['superhelical'].to_numpy()

        additional_dict['global_superhelical'] = global_superhelical
        additional_dict['local_superhelical'] = local_superhelical
        additional_dict['binding'] = binding_event
        additional_dict['unbinding'] = unbinding_event

        # Collect measurements
        # ------------------------------------------------------------------------------
        # For the enzymes.
        # ------------------------------------------------------------------------------
        kde_dict = {}
        for name in ['RNAP', 'topoI', 'gyrase']:

            # Unpack
            enzymes_df = pool_result['enzymes_df']

            # Filter
            mask = enzymes_df['name'] == name

            # Position
            # ------------------------------------------------------------------------------
            position = enzymes_df[mask]['position'].to_numpy()

            # number of bins for histogram - which we will not plot
            my_environmental = \
                [environmental for environmental in my_circuit.environmental_list if environmental.name == name][
                    0]
            if name == 'RNAP':
                size = abs(
                    my_environmental.site_list[0].start - my_environmental.size - my_environmental.site_list[0].end)
                # nbins = int(size / my_environmental.size)  # because genes are smaller
                nbins = 90  # Maybe check the one up
                x = x_gene
            else:
                nbins = int(
                    1.5*my_circuit.size / my_environmental.size)  # number of bins. - note that this only applies for topos

                x = x_system

            # Calculate KDE
            kde_x, kde_y = calculate_KDE(data=position, nbins=nbins, scaled=True)

            kde_y = get_interpolated_kde(kde_x, kde_y, x)
            kde_x = x

            # Add it to the kde_dict
            kde_dict[name] = {'kde_y': kde_y, 'kde_x': kde_x}

        # And to the additional results dict
        additional_dict['KDE'] = kde_dict

    return prod_rate, additional_dict
