import multiprocessing
from TORCphysics import parallelization_tools as pt


# Parallel stuff
n_simulations = 48  # 48 #96 #120

# Let's initialize circuit
circuit_filename = 'circuit2.csv'
sites_filename = 'sites2.csv'
enzymes_filename = 'enzymes2.csv'
environment_filename = 'environment.csv'
output_prefix = ''
frames = 2000
series = True
continuation = False
dt = 0.25

#  Run sims in parallel
Items = []
for simulation_number in range(n_simulations):
    Item = {'simulation_number': simulation_number, 'circuit_filename': circuit_filename,
            'sites_filename': sites_filename, 'enzymes_filename': enzymes_filename,
            'environment_filename': environment_filename, 'output_prefix': output_prefix,
            'frames': frames, 'series': series,
            'continuation': continuation, 'dt': dt}

    Items.append(Item)

# Create a multiprocessing pool
pool = multiprocessing.Pool()
pool_results = pool.map(pt.run_single_simulation, Items)

# colors_dict = {'tetA': '#d8b41a', 'CDS': 'silver', 'mKalama1': '#0051ff', 'Raspberry': '#e30000'}
# output = 'animation'
# out_format = '.gif'

# vs.create_animation_linear(my_circuit, my_circuit.sites_df, my_circuit.enzymes_df, output, out_format,
#                           site_type='gene', site_colours=colors_dict)
