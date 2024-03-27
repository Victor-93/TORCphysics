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


# global_dict = Dict with global simulation conditions
# variations_list = A list with a list of variations to implement to the enzymes, environmentals or sites.
# initial_substrates = A list with DNA DNA concentration
# exp_superhelicals = list with values of experimental superhelical densities
# n_simulations = how many simulations to launch.
def run_objective_function(global_dict, variations_list, exp_superhelical, n_simulations):
    # Let's run experiments for the substrate concentrations
    my_objective = 0.0
    simulation_superhelicals = []
    global_dict['DNA_concentration'] = 1  # DNA concentration - Doesn't affect

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
