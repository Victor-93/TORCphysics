#from TORCphysics.Experiments.Topo_stochastic.hyperopt_calibration.calibration_tools.calibration_tools \
#    import run_single_stochastic
from TORCphysics import Circuit
from hyperopt import tpe, hp, fmin
import numpy as np
import concurrent.futures

# import sys

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# The idea is to calibrate topo_k_cat using the hyperopt library

# ----------------------------------------------------------------------------------------------------------------------
# INITIAL CONDITIONS
# ----------------------------------------------------------------------------------------------------------------------
# Optimization conditions
num_workers = 12
chunksize = 8
nsim = 96  # in total
ntests = 500  # number of tests for parametrization
coutput = 'gyra_calibration.txt'
k_cat_min = -20.0  # Ranges to vary k_cat
k_cat_max = -5.0
alpha_min = 0.01
alpha_max = 5.0
width_min = 0.001
width_max = 5.
threshold_min = 0.001
threshold_max = 5.

k_off = 0.5
k_on = 0.005
concentration = 0.25

my_vars = ['k_cat', 'alpha', 'width', 'threshold']

# Some initial conditions
topo_k_on_0 = .005  # 0.0075
topo_k_cat_0 = -20.0
topo_concentration_0 = 0.0
gyra_k_on_0 = .005  # 0.0075
gyra_k_cat_0 = -20
gyra_concentration_0 = 0.25
gyra_k_off_0 = 0.5
topo_k_off_0 = 0.5

# Circuit conditions
circuit_filename_0 = 'circuit_3000bp_positive.csv'
sites_filename_0 = 'sites.csv'
enzymes_filename_0 = 'enzymes.csv'
environment_continuum_filename_0 = 'environment_continuum.csv'
environment_stochastic_filename_0 = 'no_bridge/single_run/environment_stochastic.csv'
tm = 'stochastic'
output_prefix = 'output'
frames = 2000
series = True
continuation = False
mm = 'uniform'
dt = .5

#BORRA ESTO
for ns in range(n_simulations):
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt, tm, mm)

    my_circuit.name = my_circuit.name + '_' + str(ns)
    my_circuit.sites_dict_list[0]['name'] = my_circuit.name
    my_circuit.log.name = my_circuit.name
    my_circuit.print_general_information()

    Lac1 = Enzyme(e_type=my_circuit.environmental_list[-1].enzyme_type,
                  name=my_circuit.environmental_list[-1].name, site=my_circuit.site_list[5],
                  position=my_circuit.site_list[5].start,
                  size=my_circuit.environmental_list[-1].size, k_cat=0.0, k_off=0.0, twist=0.0, superhelical=0.0)

    Lac2 = Enzyme(e_type=my_circuit.environmental_list[-1].enzyme_type,
                  name=my_circuit.environmental_list[-1].name, site=my_circuit.site_list[7],
                  position=my_circuit.site_list[7].start + 500,
                  size=my_circuit.environmental_list[-1].size, k_cat=0.0, k_off=0.0, twist=0.0, superhelical=0.0)

#    for site in my_circuit.site_list:
#        site.k_min = 0.01
    # This is similar to the Run function... but the idea is that we will control when the bridge is formed
    for frame in range(1, frames + 1):
        # print(frame)
        my_circuit.frame = frame
        my_circuit.time = frame * dt
        if my_circuit.series:
            my_circuit.append_sites_to_dict_step1()

        # Apply binding model and get list of new enzymes
        new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, dt, my_circuit.rng)
        if frame == bridge_time:  # Manually add the lacs
            new_enzyme_list.append(Lac1)
            new_enzyme_list.append(Lac2)
        my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

        # EFFECT
        # --------------------------------------------------------------
        effects_list = em.effect_model(my_circuit.enzyme_list, my_circuit.environmental_list, dt,
                                       my_circuit.topoisomerase_model, my_circuit.mechanical_model)
        my_circuit.apply_effects(effects_list)

        # UNBINDING
        drop_list_index, drop_list_enzyme = bm.unbinding_model(my_circuit.enzyme_list, my_circuit.dt, my_circuit.rng)
        my_circuit.drop_enzymes(drop_list_index)
        my_circuit.add_to_environment(drop_list_enzyme)

        # UPDATE GLOBALS
        my_circuit.update_global_twist()
        my_circuit.update_global_superhelical()

        # Add to series df if the series option was selected (default=True)
        if series:
            my_circuit.append_enzymes_to_dict()
            my_circuit.append_sites_to_dict_step2(new_enzyme_list, drop_list_enzyme)

    # Output the dataframes: (series)
    if series:
        my_circuit.enzymes_df = pd.DataFrame.from_dict(my_circuit.enzymes_dict_list)
        my_circuit.enzymes_df.to_csv(my_circuit.name + '_enzymes_df.csv', index=False, sep=',')
        my_circuit.sites_df = pd.DataFrame.from_dict(my_circuit.sites_dict_list)
        my_circuit.sites_df.to_csv(my_circuit.name + '_sites_df.csv', index=False, sep=',')
        my_circuit.environmental_df = pd.DataFrame.from_dict(my_circuit.environmental_dict_list)
        my_circuit.environmental_df.to_csv(my_circuit.name + '_environment_df.csv', index=False, sep=',')

    # Output the log of events
    my_circuit.log.final_twist = my_circuit.twist
    my_circuit.log.final_superhelical = my_circuit.superhelical
    my_circuit.log.log_out()

    # Output csvs
    my_circuit.enzyme_list_to_df().to_csv(my_circuit.name + '_enzymes_' + my_circuit.output_prefix + '.csv', index=False, sep=',')
    my_circuit.site_list_to_df().to_csv(my_circuit.name + '_sites_' + my_circuit.output_prefix + '.csv', index=False, sep=',')
    my_circuit.environmental_list_to_df().to_csv(my_circuit.name + '_environment_' + my_circuit.output_prefix + '.csv', index=False, sep=',')

    # And create animation
    #vs.create_animation_linear(my_circuit, my_circuit.sites_df, my_circuit.enzymes_df, my_circuit.frames,
    #                           output=my_circuit.name, out_format='.gif')


# ----------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------
def set_items(topo_kcat, topo_kon, topo_koff, topo_alpha, topo_concentration, topo_width, topo_threshold,
              gyra_kcat, gyra_kon, gyra_koff, gyra_alpha, gyra_concentration, gyra_width, gyra_threshold):
    items = []
    # for n in range(num_workers):
    for n in range(nsim):
        item = {
            'topo_concentration': topo_concentration,
            'topo_k_on': topo_kon,
            'topo_k_cat': topo_kcat,
            'topo_k_off': topo_koff,
            'topo_alpha': topo_alpha,
            'topo_width': topo_width,
            'topo_threshold': topo_threshold,
            'gyra_concentration': gyra_concentration,
            'gyra_k_on': gyra_kon,
            'gyra_k_off': gyra_koff,
            'gyra_k_cat': gyra_kcat,
            'gyra_alpha': gyra_alpha,
            'gyra_width': gyra_width,
            'gyra_threshold': gyra_threshold,
            'circuit_filename': circuit_filename_0,
            'sites_filename': sites_filename_0,
            'enzymes_filename': enzymes_filename_0,
            'environment_filename': environment_stochastic_filename_0,
            'output_prefix': output_prefix,
            'frames': frames,
            'series': series,
            'continuation': continuation,
            'dt': dt,
            'tm': 'stochastic',
            'mm': mm
        }
        items.append(item)
    return items


def objective_ProcessPoolExecutor(params):
    my_kcat = params[my_vars[0]]
    my_alpha = params[my_vars[1]]
    my_width = params[my_vars[2]]
    my_threshold = params[my_vars[3]]
    # For Gyrase
    items = set_items( topo_kcat=0.0, topo_kon=k_on, topo_koff=k_off, topo_alpha=0.0,
                       topo_concentration=0.0, topo_width=0.1, topo_threshold=0.0,
                       gyra_kcat=my_kcat, gyra_kon=k_on, gyra_koff=k_off, gyra_alpha=my_alpha,
                       gyra_concentration=concentration, gyra_width=my_width, gyra_threshold=my_threshold)
    # For Topo
    #  items = set_items( topo_kcat=my_kcat, topo_kon=k_on, topo_koff=k_off, topo_alpha=my_alpha,
    #                   topo_concentration=concentration, topo_width=my_width, topo_threshold=my_threshold,
    #                   gyra_kcat=0.0, gyra_kon=k_on, gyra_koff=k_off, gyra_alpha=0.0,
    #                   gyra_concentration=0.0, gyra_width=0.1, gyra_threshold=0.0)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit tasks to executor
        results = list(executor.map(run_single_stochastic, items, chunksize=chunksize))

    my_supercoiling = np.zeros((frames + 1, nsim))
    for i, sigma in enumerate(results):
        my_supercoiling[:, i] = sigma
    # meean = np.mean(my_supercoiling, axis=1)
    # stdd  = np.std(my_supercoiling, axis=1)
    my_objective = np.sum(np.square(np.mean(my_supercoiling, axis=1) - sigma_continuum))
    return my_objective


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
    my_vars[0]: hp.uniform(my_vars[0], k_cat_min, k_cat_max),
    my_vars[1]: hp.uniform(my_vars[1], alpha_min, alpha_max),
    my_vars[2]: hp.uniform(my_vars[2], width_min, width_max),
    my_vars[3]: hp.uniform(my_vars[3], threshold_min, threshold_max)
}

best = fmin(
    fn=objective_ProcessPoolExecutor,  # Objective Function to optimize
    space=space,  # Hyperparameter's Search Space
    algo=tpe.suggest,  # Optimization algorithm (representative TPE)
    max_evals=ntests  # Number of optimization attempts
)
print(best)
print(best.values())
best_result = np.zeros(4)
best_result[0] = best[my_vars[0]]  # Extract the best parameter
best_result[1] = best[my_vars[1]]
best_result[2] = best[my_vars[2]]
best_result[3] = best[my_vars[3]]
np.savetxt(coutput, best_result)
