from TORCphysics import Circuit
import two_genes_functions
from TORCphysics import visualization as vs
import sys
# Description: Here we will manually add RNAPs to the gene start sites so we can visualize if they behave correctly

# Initial conditions
circuit_filename_linear = 'circuit_linear.csv'
circuit_filename_circular = 'circuit_circular.csv'
sites_filename_tandem = 'sites_tandem.csv'
sites_filename_convergent = 'sites_convergent.csv'
sites_filename_divergent = 'sites_divergent.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment.csv'
output_prefix = 'output'
frames = 100
series = True
continuation = False
tm = 'continuum'
mm = 'uniform'
dt = 0.5
add_time = 10  # time in which RNAPs are added

colors_dict = {'left': '#4a86e8ff', 'right': '#ff3c3c'}

# -------------------------------------------------------------------------------------------------------------------
# LINEAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# -------------------------------------------------------------------------------------------------------------------
# Let's build the circuits
circuit_lin_tan = Circuit(circuit_filename_linear, sites_filename_tandem, enzymes_filename, environment_filename,
                          output_prefix, frames, series, continuation, dt)#, tm, mm)
circuit_lin_dir = Circuit(circuit_filename_linear, sites_filename_divergent, enzymes_filename, environment_filename,
                          output_prefix, frames, series, continuation, dt)#, tm, mm)
circuit_lin_con = Circuit(circuit_filename_linear, sites_filename_convergent, enzymes_filename, environment_filename,
                          output_prefix, frames, series, continuation, dt)#, tm, mm)

circuit_lin_tan.name = circuit_lin_tan.name + '_tandem_linear'
circuit_lin_tan.log.name = circuit_lin_tan.name
circuit_lin_dir.name = circuit_lin_dir.name + '_divergent_linear'
circuit_lin_dir.log.name = circuit_lin_dir.name
circuit_lin_con.name = circuit_lin_con.name + '_convergent_linear'
circuit_lin_con.log.name = circuit_lin_con.name
# ---------------------------------------------------------
# Tandem case
# ---------------------------------------------------------
two_genes_functions.run_simulation_two_genes(circuit_lin_tan, add_time)

#two_genes_functions.create_animation_linear(circuit_lin_tan, circuit_lin_tan.sites_df, circuit_lin_tan.enzymes_df,
#                                            frames, circuit_lin_tan.name, '.gif')
# And create visualization
my_circuit = circuit_lin_tan
#vs.create_animation_linear(my_circuit, my_circuit.sites_df, my_circuit.enzymes_df,
#                           my_circuit.name, '.gif',
#                           site_type='gene', site_colours=colors_dict)
vs.create_animation_linear_artist(my_circuit, my_circuit.sites_df, my_circuit.enzymes_df, my_circuit.environmental_df,
                           my_circuit.name, 'gif',
                           site_type='gene', site_colours=colors_dict, draw_containers=False,fps=30)
# ---------------------------------------------------------
# Divergent case
# ---------------------------------------------------------
two_genes_functions.run_simulation_two_genes(circuit_lin_dir, add_time)

#two_genes_functions.create_animation_linear(circuit_lin_dir, circuit_lin_dir.sites_df, circuit_lin_dir.enzymes_df,
#                                            frames, circuit_lin_dir.name, '.gif')
my_circuit = circuit_lin_dir
vs.create_animation_linear_artist(my_circuit, my_circuit.sites_df, my_circuit.enzymes_df, my_circuit.environmental_df,
                           my_circuit.name, 'gif',
                           site_type='gene', site_colours=colors_dict, draw_containers=False,fps=30)
# ---------------------------------------------------------
# Convergent case
# ---------------------------------------------------------
two_genes_functions.run_simulation_two_genes(circuit_lin_con, add_time)

#two_genes_functions.create_animation_linear(circuit_lin_con, circuit_lin_con.sites_df, circuit_lin_con.enzymes_df,
#                                            frames, circuit_lin_con.name, '.gif')
my_circuit = circuit_lin_con
vs.create_animation_linear_artist(my_circuit, my_circuit.sites_df, my_circuit.enzymes_df, my_circuit.environmental_df,
                           my_circuit.name, 'gif',
                           site_type='gene', site_colours=colors_dict, draw_containers=False,fps=30)

sys.exit()
# -------------------------------------------------------------------------------------------------------------------
# CIRCULAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# -------------------------------------------------------------------------------------------------------------------
# Let's build the circuits
circuit_cir_tan = Circuit(circuit_filename_circular, sites_filename_tandem, enzymes_filename, environment_filename,
                          output_prefix, frames, series, continuation, dt)
circuit_cir_dir = Circuit(circuit_filename_circular, sites_filename_divergent, enzymes_filename, environment_filename,
                          output_prefix, frames, series, continuation, dt)
circuit_cir_con = Circuit(circuit_filename_circular, sites_filename_convergent, enzymes_filename, environment_filename,
                          output_prefix, frames, series, continuation, dt)

circuit_cir_tan.name = circuit_cir_tan.name + '_tandem_circular'
circuit_cir_tan.log.name = circuit_cir_tan.name
circuit_cir_dir.name = circuit_cir_dir.name + '_divergent_circular'
circuit_cir_dir.log.name = circuit_cir_dir.name
circuit_cir_con.name = circuit_cir_con.name + '_convergent_circular'
circuit_cir_con.log.name = circuit_cir_con.name
# ---------------------------------------------------------
# Tandem case
# ---------------------------------------------------------
two_genes_functions.run_simulation_two_genes(circuit_cir_tan, add_time)

two_genes_functions.create_animation_linear(circuit_cir_tan, circuit_cir_tan.sites_df, circuit_cir_tan.enzymes_df,
                                            frames, circuit_cir_tan.name, '.gif')
# ---------------------------------------------------------
# Divergent case
# ---------------------------------------------------------
two_genes_functions.run_simulation_two_genes(circuit_cir_dir, add_time)

two_genes_functions.create_animation_linear(circuit_cir_dir, circuit_cir_dir.sites_df, circuit_cir_dir.enzymes_df,
                                            frames, circuit_cir_dir.name, '.gif')

# ---------------------------------------------------------
# Convergent case
# ---------------------------------------------------------
two_genes_functions.run_simulation_two_genes(circuit_cir_con, add_time)

two_genes_functions.create_animation_linear(circuit_cir_con, circuit_cir_con.sites_df, circuit_cir_con.enzymes_df,
                                            frames, circuit_cir_con.name, '.gif')
