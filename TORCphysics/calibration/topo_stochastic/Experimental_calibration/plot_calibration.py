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
file_out = 'calibration_1'

# For the simulation
circuit_filename = 'circuit.csv'
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
topoI_environment_filename = 'topoI_environment.csv'
gyrase_environment_filename = 'gyrase_environment.csv'

# Experimental concentration of topoisomerases
# TODO: Check bien las concentrations
gyrase_concentration = 44.6
topoI_concentration = 17.0

# Simulations conditions
output_prefix = 'test0'
series = True
continuation = False

# For parallelization and calibration
n_simulations = 2
tests = 2  # 00  # number of tests for parametrization

# Molecule/model to calibrate/test
# -----------------------------------
# Topo I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'TopoIRecognition'
topoI_effect_model_name = 'TopoIUniform'
topoI_unbinding_model_name = 'PoissonUnBinding'
topoI_calibration_file = 'topoI_calibration.csv'

# Gyrase
gyrase_name = 'topoI'
gyrase_type = 'environmental'
gyrase_binding_model_name = 'TopoIRecognition'
gyrase_effect_model_name = 'TopoIUniform'
gyrase_unbinding_model_name = 'PoissonUnBinding'
gyrase_calibration_file = 'gyrase_calibration.csv'

# FIGURE
# -----------------------------------
width = 8
height = 4
lw = 3
colors = ['blue', 'red', 'green', 'black', 'grey', 'yellow', 'cyan', 'purple', 'brown', 'pink']
fig, axs = plt.subplots(2, 2, figsize=(2 * width, 2 * height), tight_layout=True)


# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the mean superhelical
# density for each substrate concentration


def topoI_objective_function(params):
    # We need to prepare the inputs.

    # Global dictionary
    # ------------------------------------------
    global_dict = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                   'enzymes_filename': enzymes_filename, 'environment_filename': topoI_environment_filename,
                   'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                   'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': initial_sigma,
                   'DNA_concentration': 0.0}

    # Variation dictionary
    # ------------------------------------------
    name = topoI_name
    object_type = topoI_type
    binding_model_name = topoI_binding_model_name
    binding_oparams = {'k_on': params['k_on'], 'width': params['width'], 'threshold': params['threshold']}
    effect_model_name = topoI_effect_model_name
    effect_oparams = {'k_cat': params['k_cat']}
    unbinding_model_name = topoI_unbinding_model_name
    unbinding_oparams = {'k_off': params['k_off']}
    concentration = topoI_concentration

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
    binding_oparams = {'k_on': params['k_on'], 'width': params['width'], 'threshold': params['threshold']}
    effect_model_name = gyrase_effect_model_name
    effect_oparams = {'k_cat': params['k_cat']}
    unbinding_model_name = gyrase_unbinding_model_name
    unbinding_oparams = {'k_off': params['k_off']}
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
initial_sigma = -.047
final_sigma = 0.0
initial_product = 0.0
initial_substrate = .7
initial_substrates = [0.35, 0.7, 1.1, 1.8, 2.1]
enzyme_concentration = topoI_concentration
K_M = 1.5
k_cat = .0023  # 0.003
v_max = k_cat * enzyme_concentration
exp_substrates = []
exp_products = []
exp_superhelicals = []
for count, initial_substrate in enumerate(initial_substrates):
    # Substrates and products
    # ----------------------------------
    ax = axs[0, 0]
    substrate, product = tct.integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                          frames=frames, dt=dt)
    ax.plot(time, product, lw=lw, color=colors[count])

    # Sigma deduction
    # ----------------------------------
    ax = axs[1, 0]
    superhelical = tct.rescale_product_to_sigma(product, initial_sigma, final_sigma)
    ax.plot(time, superhelical, lw=lw, color=colors[count])

    # Collect
    exp_substrates.append(substrate)
    exp_products.append(product)
    exp_superhelicals.append(superhelical)

# ---------------------------------------------------------------------------------------------------------------------
# Plot simulation example for TopoI
# ---------------------------------------------------------------------------------------------------------------------
params = pd.read_csv(topoI_calibration_file).to_dict()
objective, s_superhelicals = topoI_objective_function(params=params)
name = topoI_name
object_type = topoI_type
binding_model_name = topoI_binding_model_name
binding_oparams = {'k_on': params['k_on'], 'width': params['width'], 'threshold': params['threshold']}
effect_model_name = topoI_effect_model_name
effect_oparams = {'k_cat': params['k_cat']}
unbinding_model_name = topoI_unbinding_model_name
unbinding_oparams = {'k_off': params['k_off']}
concentration = topoI_concentration

topo_variation = {'name': name, 'object_type': object_type,
                  'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                  'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                  'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                  'concentration': concentration}

objective, s_superhelicals = tct.run_objective_function(global_dict=global_dict,
                                                                    variations_list=[topo_variation],
                                                                    initial_substrates=initial_substrates,
                                                                    exp_superhelicals=exp_superhelicals,
                                                                    n_simulations=n_simulations)
# And plot
ax = axs[1, 0]
for count, initial_substrate in enumerate(initial_substrates):
    superhelical = simulation_superhelicals[count]
    ax.plot(time, superhelical, '--', lw=lw, color=colors[count])


#-------------------------------------------------------------------
for ax in [axs[0, 0], axs[0, 1]]:
    ax.set_xlabel('time (s)')
    ax.set_ylabel('Relaxed DNA (nM)')
    ax.grid(True)

for ax in [axs[1, 0], axs[1, 1]]:
    ax.set_xlabel('time (s)')
    ax.set_ylabel('Superhelical response')
    ax.grid(True)

axs[0, 0].set_title('Topoisomerase I')
axs[0, 1].set_title('Gyrase')

plt.savefig(file_out+'.png')
plt.show()
