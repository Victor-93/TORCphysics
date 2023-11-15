from TORCphysics import Circuit
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em
from TORCphysics import unbinding_model as ubm
import numpy as np


# TODO: Maybe later it can accept specific conditions
# Run simple simulation. The idea is that an external file executes this one. The external file should handle the
# parallelization process. This file is just in charge of sorting out the simulation number
def run_single_simulation(item):
    simulation_number = item['simulation_number']
    circuit_filename = item['circuit_filename']
    sites_filename = item['sites_filename']
    enzymes_filename = item['enzymes_filename']
    environment_filename = item['environment_filename']
    output_prefix = item['output_prefix']
    frames = item['frames']
    series = item['series']
    continuation = item['continuation']
    dt = item['dt']
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt)
    my_circuit.name = my_circuit.name + '_' + str(simulation_number)
    my_circuit.sites_dict_list[0]['name'] = my_circuit.name
    my_circuit.log.name = my_circuit.name
    my_circuit.run()
    return


# Helps to set the items
def set_item(circuit_filename, sites_filename, enzymes_filename, environment_filename, output_prefix, frames, series,
             continuation, dt, simulation_number):
    item = {
        'circuit_filename': circuit_filename,
        'sites_filename': sites_filename,
        'enzymes_filename': enzymes_filename,
        'environment_filename': environment_filename,
        'output_prefix': output_prefix,
        'frames': frames,
        'series': series,
        'continuation': continuation,
        'dt': dt,
        'simulation_number': simulation_number
    }
    return item


# This next functions are used for calibrating the stochastic topoisomerase model
# ----------------------------------------------------------------------------------------------------------------------
def set_items_topo_calibration(circuit_filename, sites_filename, enzymes_filename, environment_filename, output_prefix,
                               frames, series, continuation, dt, n_simulations, initial_supercoiling,
                               list_names, list_k_cat, list_k_on, list_k_off, list_width, list_threshold,
                               list_concentration, DNA_concentration):
    items = []
    for simulation_number in range(n_simulations):
        item = set_item_topo_calibration(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                                         output_prefix, frames, series, continuation, dt, simulation_number,
                                         initial_supercoiling, list_names, list_k_cat, list_k_on, list_k_off,
                                         list_width, list_threshold, list_concentration, DNA_concentration)
        items.append(item)
    return items


# Helps to set the items for calibration.
# Because we want to test different conditions for different type of enzymes, list_names contains a list of enzyme
# names that we want to vary its parameters. According to these names, the parameters on the list_* will be assigned.
def set_item_topo_calibration(circuit_filename, sites_filename, enzymes_filename, environment_filename, output_prefix,
                              frames, series, continuation, dt, simulation_number, initial_supercoiling,
                              list_names, list_k_cat, list_k_on, list_k_off, list_width, list_threshold,
                              list_concentration, DNA_concentration):
    item = {
        'circuit_filename': circuit_filename,
        'sites_filename': sites_filename,
        'enzymes_filename': enzymes_filename,
        'environment_filename': environment_filename,
        'output_prefix': output_prefix,
        'frames': frames,
        'series': series,
        'continuation': continuation,
        'dt': dt,
        'simulation_number': simulation_number,
        'initial_supercoiling': initial_supercoiling,
        'list_names': list_names,
        'list_k_cat': list_k_cat,
        'list_k_on': list_k_on,
        'list_k_off': list_k_off,
        'list_width': list_width,
        'list_threshold': list_threshold,
        'list_concentration': list_concentration,
        'DNA_concentration': DNA_concentration
    }
    return item


# This function runs another a case/system/circuit, given a set of conditions 'items'.
# def run_simulations_parallel(items, my_funciton):
#    # Create a multiprocessing pool
#    pool = multiprocessing.Pool()
#    pool_results = pool.map(my_function, items)


def run_single_simulation_topo_calibration(item):
    simulation_number = item['simulation_number']
    circuit_filename = item['circuit_filename']
    sites_filename = item['sites_filename']
    enzymes_filename = item['enzymes_filename']
    environment_filename = item['environment_filename']
    output_prefix = item['output_prefix']
    frames = item['frames']
    series = item['series']
    continuation = item['continuation']
    dt = item['dt']
    initial_supercoiling = item['initial_supercoiling']
    list_names = item['list_names']
    list_k_cat = item['list_k_cat']
    list_k_on = item['list_k_on']
    list_k_off = item['list_k_off']
    list_width = item['list_width']
    list_threshold = item['list_threshold']
    list_concentration = item['list_concentration']
    DNA_concentration = item['DNA_concentration']
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt)

    my_circuit.name = my_circuit.name + '_' + str(simulation_number)
    my_circuit.sites_dict_list[0]['name'] = my_circuit.name
    my_circuit.log.name = my_circuit.name
    #    my_circuit.superhelical = initial_supercoiling
    for enzyme in my_circuit.enzyme_list:
        enzyme.superhelical = initial_supercoiling
    my_circuit.update_twist()
    my_circuit.update_supercoiling()
    my_circuit.update_global_twist()
    my_circuit.update_global_superhelical()
    # Change topoisomerase parametrization
    for count, name in enumerate(list_names):
        for environmental in my_circuit.environmental_list:
            if environmental.name == name:
                environmental.concentration = list_concentration[count] * DNA_concentration  # Here, we are multiplying
                # [E] * [S], so the rate would be something like k = k_on * [E] * [S] * f(sigma),
                # where [E] is the enzyme concentration, [S] the substrate concentration which in this case is the
                # DNA, and f(sigma) the recognition curve.
                environmental.k_on = list_k_on[count]
                environmental.k_off = list_k_off[count]
                environmental.k_cat = list_k_cat[count]
                environmental.oparams = {'width': list_width[count], 'threshold': list_threshold[count]}

    supercoiling = my_circuit.run_return_global_supercoiling()
    return supercoiling


def single_simulation_calibration_w_supercoiling(item):
    # Item = {'global_conditions': global_dict, 'variations': variations_list}

    # Retrieve
    global_dict = item['global_conditions']
    variations_list = item['variations']

    my_circuit = Circuit(circuit_filename=global_dict['circuit_filename'], sites_filename=global_dict['sites_filename'],
                         enzymes_filename=global_dict['enzymes_filename'],
                         environment_filename=global_dict['environment_filename'],
                         output_prefix=global_dict['output_prefix'], frames=global_dict['frames'],
                         series=global_dict['series'], continuation=global_dict['continuation'],
                         dt=global_dict['dt'])

    my_circuit.name = my_circuit.name + '_' + str(global_dict['n_simulations'])
    my_circuit.sites_dict_list[0]['name'] = my_circuit.name
    my_circuit.log.name = my_circuit.name

    # And change the seed
    my_circuit.seed = my_circuit.seed + global_dict['n_simulations']  # random.randrange(sys.maxsize)
    my_circuit.rng = np.random.default_rng(my_circuit.seed)

    # Let's fix first initial supercoiling density and update all relevant parameters
    for enzyme in my_circuit.enzyme_list:
        enzyme.superhelical = global_dict['initial_sigma']
    my_circuit.update_twist()
    my_circuit.update_supercoiling()
    my_circuit.update_global_twist()
    my_circuit.update_global_superhelical()

    # And let's apply some variations.
    # Filter object to apply variations.
    # -----------------------
    for variation in variations_list:
        # We need to filter and find our my_object, which is the name of the molecule/site that we will apply the
        # variations
        if variation['object_type'] == 'enzyme':
            for enzyme in my_circuit.enzyme_list:
                if enzyme.name == variation['name']:
                    my_object = enzyme
        elif variation['object_type'] == 'environment' or variation['object_type'] == 'environmental':
            for environmental in my_circuit.environmental_list:
                if environmental.name == variation['name']:
                    my_object = environmental
                    # And let's modify concentration
                    my_object.concentration = variation['concentration'] * global_dict['DNA_concentration']
                    # Here, we are multiplying
                    # [E] * [S], so the rate would be something like k = k_on * [E] * [S] * f(sigma),
                    # where [E] is the enzyme concentration, [S] the substrate concentration which in this case is the
                    # DNA, and f(sigma) the recognition curve.

        elif variation['object_type'] == 'site':
            for site in my_circuit.site_list:
                if site.name == variation['name']:
                    my_object = site
        else:
            raise ValueError('Error, object_type not recognised')

        # Apply model variations
        # Models
        # -----------------------
        # Binding Model
        if variation['binding_model_name'] is not None:
            my_object.binding_model = bm.assign_binding_model(model_name=variation['binding_model_name'],
                                                              **variation['binding_oparams'])
        # Effect Model
        if variation['effect_model_name'] is not None:
            my_object.effect_model = em.assign_effect_model(model_name=variation['effect_model_name'],
                                                            **variation['effect_oparams'])
        # Unbinding Model
        if variation['unbinding_model_name'] is not None:
            my_object.unbinding_model = ubm.assign_unbinding_model(model_name=variation['unbinding_model_name'],
                                                                   **variation['unbinding_oparams'])

    # TODO: For some reason, the models are not passing through
    # And in case our variations are for environmentals that recognise and bind bare DNA,
    # let's define the DNA bare binding sites
    # Define bare DNA binding sites for bare DNA binding enzymes
    my_circuit.define_bare_DNA_binding_sites()
    # Sort list of enzymes and sites by position/start
    my_circuit.sort_lists()

    # Finally, run simulation
    supercoiling = my_circuit.run_return_global_supercoiling()
    return supercoiling
