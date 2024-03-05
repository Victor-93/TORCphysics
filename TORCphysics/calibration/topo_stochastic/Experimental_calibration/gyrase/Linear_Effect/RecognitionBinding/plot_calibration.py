import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import topo_calibration_tools as tct

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
topoI_environment_filename = 'topoI_environment.csv'
gyrase_environment_filename = 'gyrase_environment.csv'

# Experimental concentration of topoisomerases
# TODO: Check bien las concentrations - Creo que estan bien
gyrase_concentration = 44.6
topoI_concentration = 17.0
# Relaxed is the superhelical density when the DNA is relaxed
relaxed_DNA = 0.01
final_DNA = -0.1  # This is after the enzyme effect

# Simulations conditions
output_prefix = 'test0'
series = True
continuation = False

# For parallelization and calibration
n_simulations = 60  # 48  # 120
# Molecule/model to calibrate/test
# -----------------------------------
# Gyrase
gyrase_name = 'gyrase'
gyrase_type = 'environmental'
gyrase_binding_model_name = 'GyraseRecognition'
gyrase_effect_model_name = 'GyraseLinear'
gyrase_unbinding_model_name = 'PoissonUnBinding'
gyrase_calibration_file = 'gyrase_calibration.csv'

# FIGURE
# -----------------------------------
width = 8
height = 4
lw = 3
experiment_color = 'blue'
model_color = 'red'
fig, axs = plt.subplots(2, figsize=(width, 2 * height), tight_layout=True)


# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the mean superhelical
# density for each substrate concentration

def gyrase_objective_function(params):
    # We need to prepare the inputs.

    # Global dictionary
    # ------------------------------------------
    global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                   'enzymes_filename': enzymes_filename, 'environment_filename': gyrase_environment_filename,
                   'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                   'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': initial_sigma,
                   'DNA_concentration': 0.0}

    # Variation dictionary
    # ------------------------------------------
    name = gyrase_name
    object_type = gyrase_type
    binding_model_name = gyrase_binding_model_name
    binding_oparams = {'k_on': float(params['k_on'][0]), 'width': float(params['width'][0]),
                       'threshold': float(params['threshold'][0])}
    effect_model_name = gyrase_effect_model_name
    effect_oparams = {'k_cat': float(params['k_cat'][0]), 'sigma0': float(params['sigma0'][0])}
    unbinding_model_name = gyrase_unbinding_model_name
    unbinding_oparams = {'k_off': float(params['k_off'][0])}
    concentration = gyrase_concentration

    topo_variation = {'name': name, 'object_type': object_type,
                      'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                      'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                      'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                      'concentration': concentration}

    my_objective, simulation_superhelicals = tct.run_objective_function(global_dict=global_dict,
                                                                        variations_list=[topo_variation],
                                                                        initial_substrates=initial_substrates,
                                                                        exp_superhelicals=exp_superhelicals,
                                                                        n_simulations=n_simulations)
    return my_objective, simulation_superhelicals


# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------

# ==================================================================================================================
# ==================================================================================================================
# Build experimental curve for TOPO I
# ==================================================================================================================
# ==================================================================================================================
# Kinetics: SDNA + TopoI -> SDNA-TopoI -> RDNA + TopoI
# Product = Fluorescent or Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs
# Kinetics: SDNA + TopoI -> SDNA-TopoI -> RDNA + TopoI
# Product = Fluorescent or Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs
initial_sigma = final_DNA  # Is actually the other way around, but there's an error somewhere but I'm lazy to find it
final_sigma = relaxed_DNA  # 0.0
initial_product = 4.0
initial_substrate = .75
initial_substrates = [0.75]
enzyme_concentration = gyrase_concentration
K_M = 2.7
k_cat = .0011
v_max = k_cat * enzyme_concentration
exp_substrates = []
exp_products = []
exp_superhelicals = []
for count, initial_substrate in enumerate(initial_substrates):
    # Substrates and products
    # ----------------------------------
    ax = axs[0]
    substrate, product = tct.integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                          frames=frames, dt=dt)
    ax.plot(time, product, lw=lw, color=experiment_color, label='experiment')

    # Sigma deduction
    # ----------------------------------
    ax = axs[1]
    superhelical = tct.rescale_product_to_sigma(substrate, initial_sigma, final_sigma)
    ax.plot(time, superhelical, lw=lw, color=experiment_color, label='experiment')

    # Collect results
    exp_substrates.append(substrate)
    exp_products.append(product)
    exp_superhelicals.append(superhelical)

# ---------------------------------------------------------------------------------------------------------------------
# Plot simulation example for TopoI
# ---------------------------------------------------------------------------------------------------------------------
initial_sigma = final_sigma  # A trick
params = pd.read_csv(gyrase_calibration_file).to_dict()
objective, simulation_superhelicals = gyrase_objective_function(params=params)

# And plot
ax = axs[1]
for count, initial_substrate in enumerate(initial_substrates):
    superhelical = simulation_superhelicals[count]
    ax.plot(time, superhelical, '--', lw=lw, color=model_color, label=r'sim $\epsilon={}$'.format(objective))

# -------------------------------------------------------------------
ax = axs[0]
ax.set_xlabel('time (s)')
ax.set_ylabel('Relaxed DNA (nM)')
ax.grid(True)
ax.legend(loc='best')

ax = axs[1]
ax.set_xlabel('time (s)')
ax.set_ylabel('Superhelical response')
ax.grid(True)
ax.legend(loc='best')

axs[0].set_title('Gyrase: Linear + Recognition')

plt.savefig(file_out + '.png')
plt.savefig(file_out + '.pdf')
plt.show()
