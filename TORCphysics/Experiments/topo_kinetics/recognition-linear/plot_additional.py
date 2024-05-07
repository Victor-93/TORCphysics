import numpy as np
import matplotlib.pyplot as plt
from TORCphysics import Circuit
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em
from TORCphysics import params as papa

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This is just a test to reproduce the global supercoiling response curves from the paper:
# Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching

# TODO: Check the overflows in the binding and how much do they affect...

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 0.25
initial_time = 0
final_time = 600
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'calibration'

# For the simulation
circuit_filename = 'circuit.csv'
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
sigma_0_topo = -0.075  # Approximately -20 supercoils according the paper
sigma_0_gyrase = 0.0  # We suppose this one.
sigma_f_gyrase = -0.1  # We also assume this one, which is the maximum at which gyrase acts.
# At this value the torque is too strong.

output_prefix = 'test0'
series = True
continuation = False

# For parallelization and calibration
n_simulations = 4#84  # 60 #48 #120

# params_file
params_file = 'calibration.csv'

# Models to calibrate to calibrate
# -----------------------------------
# Topoisomerase I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'TopoIRecognition'
topoI_effect_model_name = 'TopoILinear'
topoI_unbinding_model_name = 'PoissonUnBinding'
topoI_params = 'calibration_topoI.csv'

# Gyrase
gyrase_name = 'gyrase'
gyrase_type = 'environmental'
gyrase_binding_model_name = 'GyraseRecognition'
gyrase_effect_model_name = 'GyraseLinear'
gyrase_unbinding_model_name = 'PoissonUnBinding'
gyrase_params = 'calibration_gyrase.csv'

# -----------------------------------
# FIGURE Params
# -----------------------------------
width = 8
height = 4
lw = 3
font_size = 12
title_size = 16
experiment_color = 'blue'
model_color = 'red'

# ------------------------------------------------
# Plot additional stuff: Recognition curves and effect

gyrase_color = 'blue'
topoI_color = 'red'
sigma = np.arange(-.2, .2, 0.01)

# Initialize circuit
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

# Load models because the environment file doesn't have any
# 0 = gyrase
# 1 = topoI


fig, axs = plt.subplots(2, 1, figsize=(width, 2 * height), tight_layout=True)#, sharex=True)


# Responses
# -------------------------
gyrase_bind_rate = gyrase_concentration * bm.GyraseRecognition(filename=gyrase_params).rate_modulation(superhelical=sigma)
topoI_bind_rate = topoI_concentration * bm.TopoIRecognition(filename=topoI_params).rate_modulation(superhelical=sigma)
axs[0].plot(sigma, gyrase_bind_rate, color=gyrase_color, lw=lw, label='Gyrase')
axs[0].plot(sigma, topoI_bind_rate, color=topoI_color, lw=lw, label='Topo I')

axs[0].grid(True)
axs[0].legend(loc='best')
axs[0].set_xlabel('Superhelical density')
axs[0].set_ylabel(r'Binding rate ($s^{-1}$)')

# Effects
# -------------------------
gyrase_ef = em.GyraseLinear(filename=gyrase_params).twist_added(superhelical=sigma, dt=1.0)/papa.w0
topoI_ef = em.TopoILinear(filename=topoI_params).twist_added(superhelical=sigma, dt=1.0)/papa.w0
axs[1].plot(sigma, gyrase_ef, color=gyrase_color, lw=lw)
axs[1].plot(sigma, topoI_ef, color=topoI_color, lw=lw)


axs[1].grid(True)
axs[1].set_xlabel('Superhelical density')
axs[1].set_ylabel(r'Twist added per second (bp)')

plt.savefig('additional.png')
plt.savefig('additional.pdf')
plt.show()
#my_circuit.environmental_list[3].effect_model = LacIPoissonBridging()
#my_circuit.environmental_list[3].unbinding_model = LacISimpleUnBinding()

#topoisomerase