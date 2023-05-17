from TORCphysics import Circuit
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
from hyperopt import tpe, hp, fmin
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# The idea is to calibrate topo_k_cat using the hyperopt library

# ----------------------------------------------------------------------------------------------------------------------
# INITIAL CONDITIONS
# ----------------------------------------------------------------------------------------------------------------------

# Optimization conditions
nsim = 100 # number of simulations per test
ntests = 100 #  number of tests for parametrization
# coutput = 'topoI_k_cat.txt'
coutput = 'topoI_kcat-alpha_kon.txt'
k_cat_min = 5.0  # Ranges to vary k_cat
k_cat_max = 20.0
my_var = 'topo_k_cat'
# my_vars = ['topo_k_cat', 'alpha']
my_vars = ['topo_k_cat', 'alpha', 'k_on']
alpha_min = 0.5
alpha_max = 3.0
k_on_min = 0.001
k_on_max = 0.01

# Some initial conditions
topo_k_on_0 = .005  # 0.0075
topo_k_cat_0 = 10.0
topo_concentration_0 = .25
gyra_k_on_0 = .005  # 0.0075
gyra_k_cat_0 = -20
gyra_concentration_0 = 0.0

# Circuit conditions
circuit_filename_0 = 'circuit_3000bp_negative.csv'
sites_filename_0 = 'sites.csv'
enzymes_filename_0 = 'enzymes.csv'
environment_continuum_filename_0 = 'environment_continuum.csv'
environment_stochastic_filename_0 = 'environment_stochastic.csv'
tm = 'stochastic'
output_prefix = 'output'
frames = 2000
series = True
continuation = False
mm = 'uniform'
dt = .5


# ----------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------

def objective_topoI_k_cat(params):
    my_topo_k_cat = params[my_var]
    supercoiling = run_many_stochastic(topo_concentration_0, topo_k_on_0, my_topo_k_cat,
                                       gyra_concentration_0, gyra_k_on_0, gyra_k_cat_0,
                                       circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                       environment_stochastic_filename_0, output_prefix, frames, series, continuation,
                                       dt, 'stochastic', mm)
    my_objective = np.sum(np.square(np.mean(supercoiling, axis=1) - sigma_continuum))
    return my_objective


def objective_topoI_k_cat_alpha(params):
    my_topo_k_cat = params[my_vars[0]]
    my_topo_alpha = params[my_vars[1]]
    supercoiling = run_many_stochastic(my_topo_alpha * topo_concentration_0, topo_k_on_0, my_topo_k_cat,
                                       gyra_concentration_0, gyra_k_on_0, gyra_k_cat_0,
                                       circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                       environment_stochastic_filename_0, output_prefix, frames, series, continuation,
                                       dt, 'stochastic', mm)
    my_objective = np.sum(np.square(np.mean(supercoiling, axis=1) - sigma_continuum))
    return my_objective


def objective_topoI_k_cat_alpha_k_on(params):
    my_topo_k_cat = params[my_vars[0]]
    my_topo_alpha = params[my_vars[1]]
    my_topo_k_on = params[my_vars[2]]
    supercoiling = run_many_stochastic(my_topo_alpha * topo_concentration_0, my_topo_k_on, my_topo_k_cat,
                                       gyra_concentration_0, gyra_k_on_0, gyra_k_cat_0,
                                       circuit_filename_0, sites_filename_0, enzymes_filename_0,
                                       environment_stochastic_filename_0, output_prefix, frames, series, continuation,
                                       dt, 'stochastic', mm)
    my_objective = np.sum(np.square(np.mean(supercoiling, axis=1) - sigma_continuum))
    return my_objective


def run_many_stochastic(topo_concentration, topo_k_on, topo_k_cat, gyra_concentration, gyra_k_on, gyra_k_cat,
                        circuit_filename, sites_filename, enzymes_filename, environment_filename,
                        output_prefix, nframes, series, continuation, dt, tm, mm):
    supercoiling = np.zeros((nframes + 1, nsim))
    for i in range(nsim):
        supercoiling[:, i] = run_stochastic_sim(topo_concentration, topo_k_on, topo_k_cat,
                                                gyra_concentration, gyra_k_on, gyra_k_cat,
                                                circuit_filename, sites_filename, enzymes_filename,
                                                environment_filename, output_prefix, nframes, series, continuation,
                                                dt, tm, mm)
    return supercoiling


def run_stochastic_sim(topo_concentration, topo_k_on, topo_k_cat, gyra_concentration, gyra_k_on, gyra_k_cat,
                       circuit_filename, sites_filename, enzymes_filename, environment_filename,
                       output_prefix, nframes, series, continuation, dt, tm, mm):
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, nframes, series, continuation, dt, tm, mm)
    my_circuit.environmental_list[0].k_on = topo_k_on
    my_circuit.environmental_list[0].k_cat = topo_k_cat
    my_circuit.environmental_list[0].concentration = topo_concentration
    my_circuit.environmental_list[1].k_on = gyra_k_on
    my_circuit.environmental_list[1].k_cat = gyra_k_cat
    my_circuit.environmental_list[1].concentration = gyra_concentration

    my_supercoiling = np.zeros(nframes + 1)
    my_supercoiling[0] = my_circuit.superhelical

    # run simulation
    for frame in range(1, nframes + 1):
        # BINDING
        # --------------------------------------------------------------
        new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                           my_circuit.rng)
        my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

        # EFFECT
        # --------------------------------------------------------------
        effects_list = em.effect_model(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                       my_circuit.topoisomerase_model, my_circuit.mechanical_model)
        my_circuit.apply_effects(effects_list)

        # UNBINDING
        # --------------------------------------------------------------
        drop_list_index, drop_list_enzyme = bm.unbinding_model(my_circuit.enzyme_list, my_circuit.dt,
                                                               my_circuit.rng)
        my_circuit.drop_enzymes(drop_list_index)
        my_circuit.add_to_environment(drop_list_enzyme)
        # UPDATE GLOBALS
        # --------------------------------------------------------------
        my_circuit.update_global_twist()
        my_circuit.update_global_superhelical()
        my_supercoiling[frame] = my_circuit.superhelical
    return my_supercoiling


# ----------------------------------------------------------------------------------------------------------------------
# PROCESS
# ----------------------------------------------------------------------------------------------------------------------

# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------
tm = 'continuum'
continuum_circuit = Circuit(circuit_filename_0, sites_filename_0, enzymes_filename_0, environment_continuum_filename_0,
                            output_prefix, frames, series, continuation, dt, tm, mm)
continuum_circuit.environmental_list[1].concentration = gyra_concentration_0
continuum_circuit.run()

# Get global supercoiling responses from the continuum case
mask = continuum_circuit.sites_df['type'] == 'circuit'  # This one contains global superhelical density
sigma_continuum = continuum_circuit.sites_df[mask]['superhelical'].to_numpy()

# Optimization case
# ----------------------------------------------------------------------------------------------------------------------
space = {
#    my_var: hp.uniform(my_var, k_cat_min, k_cat_max),  # THIS WAS ONLY WITH TOPO K_CAT
    my_vars[0]: hp.uniform(my_vars[0], k_cat_min, k_cat_max),
    my_vars[1]: hp.uniform(my_vars[1], alpha_min, alpha_max),
    my_vars[2]: hp.uniform(my_vars[2], k_on_min, k_on_max)
}
best = fmin(
#    fn=objective_topoI_k_cat_alpha,  # Objective Function to optimize
    fn=objective_topoI_k_cat_alpha_k_on,  # Objective Function to optimize
    space=space,  # Hyperparameter's Search Space
    algo=tpe.suggest,  # Optimization algorithm (representative TPE)
    max_evals=ntests  # Number of optimization attempts
)
print(best)
print(best.values())
best_result = np.zeros(3)
#best_result[0] = best[my_var]  # Extract the best parameter
best_result[0] = best[my_vars[0]]  # Extract the best parameter
best_result[1] = best[my_vars[1]]  # Extract the best parameter
best_result[2] = best[my_vars[2]]  # Extract the best parameter
np.savetxt(coutput, best_result)
