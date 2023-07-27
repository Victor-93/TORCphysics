from TORCphysics import Circuit
import multiprocessing


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
    tm = item['tm']
    mm = item['mm']
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, frames, series, continuation, dt, tm, mm)
    my_circuit.name = my_circuit.name + '_' + str(simulation_number)
    my_circuit.sites_dict_list[0]['name'] = my_circuit.name
    my_circuit.log.name = my_circuit.name
    my_circuit.run()
    return


# Helps to set the items
def set_item(circuit_filename, sites_filename, enzymes_filename, environment_filename, output_prefix, frames, series,
             continuation, dt, tm, mm, simulation_number):
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
        'tm': tm,
        'mm': mm,
        'simulation_number': simulation_number
    }
    return item


# This next functions are used for calibrating the stochastic topoisomerase model
# ----------------------------------------------------------------------------------------------------------------------
def set_items_topo_calibration(circuit_filename, sites_filename, enzymes_filename, environment_filename, output_prefix,
                               frames, series, continuation, dt, tm, mm, n_simulations, initial_supercoiling,
                               list_names, list_k_cat, list_k_on, list_k_off, list_width, list_threshold,
                               list_concentration, DNA_concentration):
    items = []
    for simulation_number in range(n_simulations):
        item = set_item_topo_calibration(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                                         output_prefix, frames, series, continuation, dt, tm, mm, simulation_number,
                                         initial_supercoiling, list_names, list_k_cat, list_k_on, list_k_off,
                                         list_width, list_threshold, list_concentration, DNA_concentration)
        items.append(item)
    return items


# Helps to set the items for calibration.
# Because we want to test different conditions for different type of enzymes, list_names contains a list of enzyme
# names that we want to vary its parameters. According to these names, the parameters on the list_* will be assigned.
def set_item_topo_calibration(circuit_filename, sites_filename, enzymes_filename, environment_filename, output_prefix,
                              frames, series, continuation, dt, tm, mm, simulation_number, initial_supercoiling,
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
        'tm': tm,
        'mm': mm,
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
#def run_simulations_parallel(items, my_funciton):
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
    tm = item['tm']
    mm = item['mm']
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
                         output_prefix, frames, series, continuation, dt, tm, mm)
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
