from TORCphysics import Circuit

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
