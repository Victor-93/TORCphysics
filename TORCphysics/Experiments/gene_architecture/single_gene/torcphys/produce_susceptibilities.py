from TORCphysics import Circuit
import numpy as np
import multiprocessing
import pandas as pd
from TORCphysics import parallelization_tools as pt

# TODO: Create template circuit and sites. Change these conditions per system conditions. - You can create temporal
#       circuit and sites files per system. Then create new ones and overwrite.
#       Run in parallel each of the systems.
#       For each run, collect only data of interest (susceptibilities), or you could save as binary to study dynamics.
#       Quantities worth of getting would be production rates (unbinding), global supercoiling, supercoiling at promoter,
#       positions? Could be the case. Maybe KDEs like you did in RNAPTracking (so you don't accomulate too much data).

# **********************************************************************************************************************
# Inputs/Initial conditions
# **********************************************************************************************************************

# Simulation conditions
# --------------------------------------------------------------
file_out = 'output'
dt = 0.25
initial_time = 0
final_time = 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Parallelization conditions
# --------------------------------------------------------------
n_simulations = 4

# Circuit conditions
# --------------------------------------------------------------
output_prefix = 'single_gene'  #Although we won't use it
series = True
continuation = False
enzymes_filename = None
environmental_filename = ''
# We will define sites and circuit csv as we advance

# Junier experimental data - we only need distances for now.
# --------------------------------------------------------------
weak_exp = '../junier_data/weak.csv'
medium_exp = '../junier_data/medium.csv'
strong_exp = '../junier_data/strong.csv'

# Promoter responses
# --------------------------------------------------------------
weak_response = '../promoter_data/weak.csv'
medium_response = '../promoter_data/medium.csv'
strong_response = '../promoter_data/strong.csv'


# **********************************************************************************************************************
# Functions
# **********************************************************************************************************************
# This functions need to create the csv files and write them
def make_linear_circuit_csv(filename, nbp, sigma0, name):
    info = {'name': name, 'structure': 'linear', 'size': nbp, 'twist': 0.0, 'superhelical': sigma0, 'sequence': None}
    circuit_df = pd.DataFrame([info])
    circuit_df.to_csv(filename, index=False)


# Make circuit
# Make site

# **********************************************************************************************************************
# Load data and prepare systems
# **********************************************************************************************************************

# Load exp
weak_exp = pd.read_csv(weak_exp)
medium_exp = pd.read_csv(medium_exp)
strong_exp = pd.read_csv(strong_exp)

# Extract distances
weak_distances = list(weak_exp['distance'])
medium_distances = list(medium_exp['distance'])
strong_distances = list(strong_exp['distance'])

# TODO: Test writing circuit
# TODO: The size is not the distance of the upstream barrier.
#       It should be d_upstream + gene_length + d_downstream. - Define this variables then make it
make_linear_circuit_csv('circuit_test.csv', weak_distances[0], -0.04, 'weak')
