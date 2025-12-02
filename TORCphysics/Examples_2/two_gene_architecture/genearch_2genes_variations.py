import numpy as np
from TORCphysics import Circuit, Site
import random
import pickle
import sys
from itertools import product

from TORCphysics.src import model_params_dir # This one is to load promoter parameterisations

# **********************************************************************************************************************
# Description
# **********************************************************************************************************************
# THIS IS THE THIRD EXAMPLE
# Based on the experiments of assessing the impact of genome architecture...
# This one tries to add a second gene. At the moment both genes are according the promoter_case (try mixing them).
# Here, we will try varying domain sizes, gene orientation, and promoters.

# TODO: 1) Try adding a second gene and modifying the orientation. 2) Mixing promoters
#       3) Increasing the distances between the gene, and the barriers.
#

# **********************************************************************************************************************
# Inputs/Initial conditions - At least the ones you need to change before each run
# **********************************************************************************************************************
dt=1.0
initial_time = 0
# final_time = 5400 #~1.5hrs
final_time = 500 # 500 is ok for testing.
n_simulations = 4 # Number of simulations per system

RNAP_gamma = 0.1

file_out = 'genearch_example_3'

time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Initial superhelical density
sigma0 = -0.046

# Model that describes the binding of promoters/genes
genes_binding_model = 'GaussianBinding'

# Circuit files
environment_filename = 'environment.csv'

# System variation conditions
# --------------------------------------------------------------
# The system is defined as: B__DIST1_____GENE1_____INTERLENGTH____GENE2____DIST2_____B
# where:
#   B: barriers at the end (boundaries)
#   DIST1 : Distance between left barrier and gene on the left (or gene1)
#   GENE1 : Length of gene 1 (left)
#   INTERLENGTH : Distance between genes
#   GENE2 : Length of gene 2 (right)
#   DIST3 : Distance between gene 2 and barrier on the right

# Length - We can vary the lengths of genes and distances between them and barriers
#dist1_list = [320, 1000]
dist1_list = [1000]
intergene_length_list = [200,1000]  # Only variation right now because it'd be massive
dist2_list = [1000]
gene_length_list = [900]  # We may not want to vary this one - for now, this will be the length of both genes

# Orientations - We can vary gene orientations
#gene1_orientation = [-1, 1]
gene1_orientation = [-1] # Indicates the orientation. If -1, then it points to the left.
                         # If 1, then points to the right
gene2_orientation = [1]

# Promoters - We can vary promoters as well. Let's use promoters that we already calibrated in the TORCphysics paper.
gene1_promoters = ['medium'] #['weak', 'medium']
gene2_promoters = ['strong']  #[0,1,2]

total_number_of_systems = (len(dist1_list) * len(intergene_length_list) * len(dist2_list) * len(gene_length_list) *
                           len(gene1_orientation) *len(gene2_orientation) * len(gene1_promoters) * len(gene2_promoters))
print('Total number of systems that will be simulated: ', total_number_of_systems)
print('Number of simulations per system: ', n_simulations)
print('Total number of simulations: ', total_number_of_systems*n_simulations)

# **********************************************************************************************************************
# Prepare system and run it
# **********************************************************************************************************************

system_outputs = [] # Here, we will store all the results

for dist1, intergene_length, dist2, gene_length, g1_promoter, g2_promoter in product(
        dist1_list, intergene_length_list, dist2_list, gene_length_list, gene1_promoters, gene2_promoters):

    # Build the circuit
    # ----------------------------------------------
    csize = dist1+gene_length+intergene_length+gene_length+dist2 # This because: B__DIST1_____GENE1_____INTERLENGTH____GENE2____DIST2_____B
    my_circuit = Circuit(name='example3', structure='linear', superhelicity=sigma0, size=csize,
                         environment_filename=environment_filename, frames=frames, dt=dt)
    # We just gave it the environmental csv file, but we can also define it through scripting.

    # Build Sites
    # ----------------------------------------------
    g1_start = dist1
    g1_end = g1_start + gene_length
    g1_params = model_params_dir + '/promoter_' +  g1_promoter + '.csv' # These parameters are pre-saved!

    g2_start = dist1 + gene_length +  intergene_length
    g2_end = g2_start + gene_length
    g2_params = model_params_dir + '/promoter_' + g2_promoter + '.csv'

    # Define the sites and add them to the circuit
    gene1 = Site(site_type='gene', name='left', start=g1_start, end=g1_end,
                 binding_model_name=genes_binding_model, binding_oparams_file=g1_params)
    gene2 = Site(site_type='gene', name='right', start=g2_start, end=g2_end,
                 binding_model_name=genes_binding_model, binding_oparams_file=g2_params)

    my_circuit.add_custom_Site(gene1)  # This function add the sites to the circuit.
    my_circuit.add_custom_Site(gene2)

    # Run simulations!
    # -----------------------------------------------------------------------------
    # Prepare lists that will store the outputs for each condition
    enzymes_df_list = []
    sites_df_list = []
    environment_df_list = []

    for n in range(n_simulations): # Run multiple simulations

        my_circuit.name = my_circuit.name + '_' + str(n)  # We can change the name of the circuit, so it is associated with each sim
        my_circuit.seed = my_circuit.seed + n + random.randrange(sys.maxsize)  # Let's change the seed!
        my_circuit.rng = np.random.default_rng(my_circuit.seed)

        # Run simulations but storing dataframes on memory then adding them to the lists.
        enzymes_df, sites_df, environment_df = my_circuit.run_return_dfs()  # Function run_return_dfs() returns the dataframes with the results of the simulation (it does not write CSV files).

        # Append dataframes to the lists.
        enzymes_df_list.append(enzymes_df)
        sites_df_list.append(sites_df)
        environment_df_list.append(environment_df)

    internal_results_dict = {'sites_df': sites_df_list, 'enzymes_df': enzymes_df_list, 'environment_df': environment_df_list}

    # Let's prepare the output of the current system iteration.
    # -----------------------------------------------------------------------------
    # The idea is to have a dictionary for each built system.
    # We store information of the build system like sizes, promoters, etc...
    # And store the results from the simulaitons

    info_dict = {'dist1': dist1, 'intergene_length': intergene_length, 'dist2': dist2, 'gene_length': gene_length,
                    'g1_promoter':g1_promoter, 'g2_promoter':g2_promoter}
    system_outputs.append({'system_info': info_dict, 'results': internal_results_dict})

# Save the dictionary to a file
with open(file_out + '.pkl', 'wb') as file:
    pickle.dump(system_outputs, file)

# TODO: Try plotting the results!
