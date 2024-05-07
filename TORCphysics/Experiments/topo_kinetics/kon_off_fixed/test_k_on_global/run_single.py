import sys
import numpy as np
from TORCphysics import Circuit

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This is to reproduce the global supercoiling response curves from the paper:
# Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching

# This time will deduce k_on & k_off based on the kinetic measurements and assuming that k_off = 0.5 (2 seconds),
# based on:
# Single-molecule imaging of DNA gyrase activity in living Escherichia coli

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 0.25
initial_time = 0
# Let's do it for 400s to add more weight to the curve and not the plateau
final_time = 5000
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'calibration'

# For the simulation
circuit_filename = '../../circuit.csv'
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
environment_filename = 'environment.csv'

# Concentrations in nM
DNA_concentration = 0.75
gyrase_concentration = 44.6
topoI_concentration = 17.0

# MM kinetics
K_M_topoI = 1.5
k_cat_topoI = .0023
v_max_topoI = k_cat_topoI * topoI_concentration
K_M_gyrase = 2.7
k_cat_gyrase = .0011
v_max_gyrase = k_cat_gyrase * gyrase_concentration

# Superhelical values (sigma) for each case
sigma_0_topo = -0.11  #-0.076  # Approximately -20 supercoils according the paper
sigma_0_gyrase = 0.0  # We suppose this one.
sigma_f_gyrase = -0.11  # We also assume this one, which is the maximum at which gyrase acts.
# At this value the torque is too strong.

output_prefix = 'test0'
series = True
continuation = False

# Models to calibrate to calibrate
# -----------------------------------
# Topoisomerase I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'TopoIRecognition'
topoI_effect_model_name = 'TopoILinear'
topoI_unbinding_model_name = 'PoissonUnBinding'

# Gyrase
gyrase_name = 'gyrase'
gyrase_type = 'environmental'
gyrase_binding_model_name = 'GyraseRecognition'
gyrase_effect_model_name = 'GyraseLinear'
gyrase_unbinding_model_name = 'PoissonUnBinding'

# SOME SIZES
# -----------------------------------
plasmid_length = 2757 # Make sure is the same as in circuit.csv
topoI_size = 60
gyrase_size = 160
n_grid_topoI = int(plasmid_length/topoI_size)
n_grid_gyrase = int(plasmid_length/gyrase_size)

# FIXED VALUES
# -----------------------------------
k_off_topoI = 0.5
k_off_gyrase = 0.5

k_on_topoI_exp = (k_off_topoI + k_cat_topoI) / K_M_topoI
k_on_gyrase_exp = (k_off_gyrase + k_cat_gyrase) / K_M_gyrase

# Deduce k_on based on k_off
k_on_topoI = ((k_off_topoI + k_cat_topoI) / K_M_topoI)/(n_grid_topoI*topoI_concentration)
k_on_gyrase = ((k_off_gyrase + k_cat_gyrase) / K_M_gyrase)/(n_grid_gyrase*gyrase_concentration)

print('k_on topoI:', k_on_topoI, 'k_on_exp topoI:', k_on_topoI_exp)
print('k_on_gyrase:', k_on_gyrase, 'k_on_exp gyrase:', k_on_gyrase_exp)

my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)
my_circuit.print_general_information()
my_circuit.run()

