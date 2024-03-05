import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from TORCphysics import binding_model as bm
from TORCphysics import Circuit
# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# We want to plot the results of the calibration

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 0.5
initial_time = 0
final_time = 600
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'calibration'

# For the simulation
circuit_filename = 'circuit.csv'
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
environment_filename = 'topoI_environment.csv'
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     'output', frames, True, False, dt)

# Experimental concentration of topoisomerases
# TODO: Check bien las concentrations - Creo que estan bien
gyrase_concentration = 44.6
topoI_concentration = 17.0
my_circuit.environmental_list[0].concentration = topoI_concentration
topoI_environmental = my_circuit.environmental_list[0]

# Simulations conditions
output_prefix = 'test0'
series = True
continuation = False

# For parallelization and calibration
n_simulations = 48 #120
# Molecule/model to calibrate/test
# -----------------------------------
# Topo I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'TopoIRecognition'
topoI_effect_model_name = 'TopoisomeraseLinearEffect'
#topoI_effect_model_name = 'TopoisomeraseLinearRandEffect'
topoI_unbinding_model_name = 'PoissonUnBinding'
topoI_calibration_file = 'topoI_calibration.csv'
topoI_sigma0 = 0.0

params = pd.read_csv(topoI_calibration_file).to_dict()

name = topoI_name
object_type = topoI_type
binding_model_name = topoI_binding_model_name
binding_oparams = {'k_on': float(params['k_on'][0]), 'width': float(params['width'][0]),
                   'threshold': float(params['threshold'][0])}
effect_model_name = topoI_effect_model_name
effect_oparams = {'k_cat': float(params['k_cat'][0]), 'sigma0': topoI_sigma0}
unbinding_model_name = topoI_unbinding_model_name
unbinding_oparams = {'k_off': float(params['k_off'][0])}
concentration = topoI_concentration

# FIGURE
# -----------------------------------
width = 8
height = 4
lw = 3
model_color = 'red'
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)

sigma = np.arange(-.1, .1, .001)
rate = np.zeros_like(sigma)
binding_model = bm.assign_binding_model(model_name=binding_model_name, **binding_oparams)
for i, si in enumerate(sigma):
    rate[i] = binding_model.binding_probability(topoI_environmental, si,dt)

axs.plot(sigma, rate, color=model_color, lw=2)
axs.set_xlabel('sigma')
axs.set_ylabel('probability')
axs.grid(True)
plt.savefig('probability_response.png')
