from TORCphysics import Circuit
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import numpy as np


def run_single_stochastic(item):
    topo_concentration = item['topo_concentration']
    topo_k_on = item['topo_k_on']
    topo_k_cat = item['topo_k_cat']
    topo_alpha = item['topo_alpha']
    topo_k_off = item['topo_k_off']
    topo_width = item['topo_width']
    topo_threshold = item['topo_threshold']
    gyra_concentration = item['gyra_concentration']
    gyra_k_on = item['gyra_k_on']
    gyra_k_off = item['gyra_k_off']
    gyra_k_cat = item['gyra_k_cat']
    gyra_alpha = item['gyra_alpha']
    gyra_width = item['gyra_width']
    gyra_threshold = item['gyra_threshold']
    circuit_filename = item['circuit_filename']
    sites_filename = item['sites_filename']
    enzymes_filename = item['enzymes_filename']
    environment_filename = item['environment_filename']
    output_prefix = item['output_prefix']
    nframes = item['frames']
    series = item['series']
    continuation = item['continuation']
    dt = item['dt']
    tm = item['tm']
    mm = item['mm']
    supercoiling = run_stochastic_sim(topo_concentration, topo_k_on, topo_k_cat, topo_k_off,
                                      topo_alpha, topo_width, topo_threshold,
                                      gyra_concentration, gyra_k_on, gyra_k_cat, gyra_k_off,
                                      gyra_alpha, gyra_width, gyra_threshold,
                                      circuit_filename, sites_filename, enzymes_filename,
                                      environment_filename, output_prefix, nframes, series, continuation,
                                      dt, tm, mm)
    return supercoiling


def run_many_stochastic(item, nsim=10):
    topo_concentration = item['topo_concentration']
    topo_k_on = item['topo_k_on']
    topo_k_cat = item['topo_k_cat']
    gyra_concentration = item['gyra_concentration']
    gyra_k_on = item['gyra_k_on']
    gyra_k_cat = item['gyra_k_cat']
    gyra_width = item['gyra_width']
    gyra_threshold = item['gyra_threshold']
    circuit_filename = item['circuit_filename']
    sites_filename = item['sites_filename']
    enzymes_filename = item['enzymes_filename']
    environment_filename = item['environment_filename']
    output_prefix = item['output_prefix']
    nframes = item['frames']
    series = item['series']
    continuation = item['continuation']
    dt = item['dt']
    tm = item['tm']
    mm = item['mm']
    supercoiling = np.zeros((nframes + 1, nsim))
    for i in range(nsim):
        supercoiling[:, i] = run_stochastic_sim(topo_concentration, topo_k_on, topo_k_cat,
                                                gyra_concentration, gyra_k_on, gyra_k_cat,
                                                circuit_filename, sites_filename, enzymes_filename,
                                                environment_filename, output_prefix, nframes, series, continuation,
                                                dt, tm, mm)
    return supercoiling


def run_stochastic_sim(topo_concentration, topo_k_on, topo_k_cat, topo_k_off,
                       topo_alpha, topo_width, topo_threshold,
                       gyra_concentration, gyra_k_on, gyra_k_cat, gyra_k_off,
                       gyra_alpha, gyra_width, gyra_threshold,
                       circuit_filename, sites_filename, enzymes_filename, environment_filename,
                       output_prefix, nframes, series, continuation, dt, tm, mm):

    # Initialize circuit
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, nframes, series, continuation, dt, tm, mm)

    # Change params
    my_circuit.environmental_list[0].concentration = topo_concentration
    my_circuit.environmental_list[0].k_on = topo_k_on * topo_alpha
    my_circuit.environmental_list[0].k_off = topo_k_off
    my_circuit.environmental_list[0].k_cat = topo_k_cat
    my_circuit.environmental_list[0].oparams = {'width': topo_width, 'threshold': topo_threshold}

    my_circuit.environmental_list[1].concentration = gyra_concentration
    my_circuit.environmental_list[1].k_on = gyra_k_on * gyra_alpha
    my_circuit.environmental_list[1].k_cat = gyra_k_cat
    my_circuit.environmental_list[1].k_off = gyra_k_off
    my_circuit.environmental_list[1].oparams = {'width': gyra_width, 'threshold': gyra_threshold}

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
