import numpy as np
import pandas as pd
import topo_calibration_tools as tct
import matplotlib.pyplot as plt

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
file_out = 'calibration_Poisson'

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
n_simulations = 84  # 60 #48 #120

# params_file
params_file = 'calibration.csv'

# Models to calibrate to calibrate
# -----------------------------------
# Topoisomerase I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'PoissonBinding'
topoI_effect_model_name = 'TopoILinear'
topoI_unbinding_model_name = 'PoissonUnBinding'
topoI_params = 'calibration_topoI.csv'

# Gyrase
gyrase_name = 'gyrase'
gyrase_type = 'environmental'
gyrase_binding_model_name = 'PoissonBinding'
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
xlabel_size = 14
title_size = 16
experiment_color = 'blue'
model_color = 'red'


# Optimization functions
# ----------------------------------------------------------------------------------------------------------------------
# This one runs the objective function in parallel. It returns the objective function as well as the mean superhelical
# density for each substrate concentration

def objective_function(params):
    # We need to prepare the inputs.
    # This time we have three different systems:
    # 1.- Topoisomerase I acting on supercoiled DNA
    # 2.- Gyrase acting on Relaxed DNA
    # 3.- Both enzymes acting on supercoiled/Relaxed DNA

    # Global dictionaries
    # ------------------------------------------
    global_dict_topoI = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                         'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                         'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                         'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma_0_topo,
                         'DNA_concentration': 0.0}

    global_dict_gyrase = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                          'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                          'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                          'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma_0_gyrase,
                          'DNA_concentration': 0.0}

    global_dict_both = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                        'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                        'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                        'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma_0_topo,
                        'DNA_concentration': 0.0}

    # Variation dictionaries
    # ------------------------------------------

    # Topoisomerase I
    name = topoI_name
    object_type = topoI_type
    binding_model_name = topoI_binding_model_name
    binding_oparams = {'k_on': float(params['k_on_topoI'][0])}
    effect_model_name = topoI_effect_model_name
    effect_oparams = {'k_cat': float(params['k_cat_topoI'][0])}
    unbinding_model_name = topoI_unbinding_model_name
    unbinding_oparams = {'k_off': float(params['k_off_topoI'][0])}
    concentration = topoI_concentration  # / mol_concentration  # Because this is the reference.

    topoI_variation = {'name': name, 'object_type': object_type,
                       'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                       'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                       'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                       'concentration': concentration}

    # Gyrase
    name = gyrase_name
    object_type = gyrase_type
    binding_model_name = gyrase_binding_model_name
    binding_oparams = {'k_on': float(params['k_on_gyrase'][0])}
    effect_model_name = gyrase_effect_model_name
    effect_oparams = {'k_cat': float(params['k_cat_gyrase'][0]), 'sigma0': float(params['sigma0_gyrase'][0])}
    unbinding_model_name = gyrase_unbinding_model_name
    unbinding_oparams = {'k_off': float(params['k_off_gyrase'][0])}
    concentration = gyrase_concentration  # / mol_concentration  # Because this is the reference.

    gyrase_variation = {'name': name, 'object_type': object_type,
                        'binding_model_name': binding_model_name, 'binding_oparams': binding_oparams,
                        'effect_model_name': effect_model_name, 'effect_oparams': effect_oparams,
                        'unbinding_model_name': unbinding_model_name, 'unbinding_oparams': unbinding_oparams,
                        'concentration': concentration}

    # Create lists of conditions for each system
    # ------------------------------------------

    # Global dictionaries
    global_dict_list = [global_dict_topoI, global_dict_gyrase, global_dict_both]

    # List of lists of variations
    variations_list = [[topoI_variation], [gyrase_variation], [topoI_variation, gyrase_variation]]

    # Arrays with global superhelical densities
    list_sigmas = [topoI_sigma, gyrase_sigma, both_sigma]

    # Finally, run objective function. run_objective_function will process our conditions
    # ------------------------------------------
    my_objective, simulation_superhelicals = tct.run_objective_function(global_dict_list=global_dict_list,
                                                                        variations_list=variations_list,
                                                                        exp_superhelicals=list_sigmas,
                                                                        n_simulations=n_simulations)
    return my_objective, simulation_superhelicals


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Experimental curves
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------
# Build experimental curve for TOPO I
# -----------------------------------------------------
# Kinetics: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI
# Product = Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs; which initially is the same as the DNA concentration

# Integrate MM kinetics
# ------------------------------------------
# Initially, there's no relaxed DNA, and all the supercoiled DNA concentration corresponds to the plasmid conc.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_topoI(vmax=v_max_topoI, KM=K_M_topoI,
                                                      Supercoiled_0=DNA_concentration, Relaxed_0=0.0,
                                                      frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
topoI_sigma = tct.topoI_to_sigma(Relaxed=relaxed_DNA, DNA_concentration=DNA_concentration, sigma0=sigma_0_topo)

# -----------------------------------------------------
# Build experimental curve for Gyrase
# -----------------------------------------------------
# Kinetics: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Product = Supercoiled DNA
# Substrate = Relaxed DNA; which initially is the same as the DNA concentration

# Integrate MM kinetics
# ------------------------------------------
# Initially, there's no supercoiled DNA, and all of the relaxed DNA concentration corresponds
# to the plasmid concentration.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_gyrase(vmax=v_max_gyrase, KM=K_M_gyrase,
                                                       Supercoiled_0=0.0, Relaxed_0=DNA_concentration,
                                                       frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
gyrase_sigma = tct.gyrase_to_sigma(Relaxed=relaxed_DNA, DNA_concentration=DNA_concentration,
                                   sigma0=sigma_0_gyrase, sigmaf=sigma_f_gyrase)

# -----------------------------------------------------
# Build experimental curve for system with both Topo I and Gyrase
# -----------------------------------------------------
# Kinetics Gyrase: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Kinetics Topoisomerase: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI

# Integrate MM kinetics
# ------------------------------------------

# Initially, there's no supercoiled DNA, and all of the relaxed DNA concentration corresponds
# to the plasmid concentration.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_both_T_G(vmax_topoI=v_max_topoI, vmax_gyrase=v_max_gyrase,
                                                         KM_topoI=K_M_topoI, KM_gyrase=K_M_gyrase,
                                                         Supercoiled_0=DNA_concentration, Relaxed_0=0.0,
                                                         frames=frames, dt=dt)
ratio = relaxed_DNA[-1] / DNA_concentration
sigmaf = sigma_0_topo * ratio
# Translate to superhelical density
# ------------------------------------------
both_sigma = tct.both_T_G_to_sigma(Relaxed=relaxed_DNA, Relaxed_final=relaxed_DNA[-1],
                                   sigma0=sigma_0_topo, sigmaf=sigmaf)

# ----------------------------------------------------------------------------------------------------------------------
# Theoretical curves
# ----------------------------------------------------------------------------------------------------------------------

params_dict = pd.read_csv(params_file).to_dict()
objective, sim_superhelicals = objective_function(params=params_dict)

# ----------------------------------------------------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------------------------------------------------

# Create figure
fig, axs = plt.subplots(3, 1, figsize=(width, 3 * height), tight_layout=True)

# Prepare arrays to iterate
exp_sigmas = [topoI_sigma, gyrase_sigma, both_sigma]
titles = ['Topoisomerase I', 'Gyrase', ' Topoisomerase & Gyrase']
outside_label = ['a)', 'b)', 'c)']
# Plot
# ------------------------------------------
for n in range(3):
    ax = axs[n]
    sum_err = np.sum(np.square(sim_superhelicals[n] - exp_sigmas[n]))

    if n == 0:
        ax.plot(time, exp_sigmas[n], lw=lw, color=experiment_color, label='kinetic')
        ax.plot(time, sim_superhelicals[n], lw=lw, color=model_color, label='model')
        ax.legend(loc='best', fontsize=font_size)
    else:
        ax.plot(time, exp_sigmas[n], lw=lw, color=experiment_color)
        ax.plot(time, sim_superhelicals[n], lw=lw, color=model_color)

    ax.set_xlabel('Time (s)', fontsize=xlabel_size)
    ax.set_ylabel('Superhelical density', fontsize=xlabel_size)
    ax.grid(True)
    ax.set_title(titles[n], fontsize=title_size)

    # Add errors
    # Create the label
    formatted_sum_err = "{:.3g}".format(sum_err)  # Formatting to three significant figures
    er_label = r'$\epsilon={}$'.format(formatted_sum_err)
    text_x = 0.85
    text_y = 0.5
    # Add the text to the plot
    ax.text(text_x, text_y, er_label, fontsize=font_size, transform=ax.transAxes,
            bbox=dict(facecolor='gray', alpha=0.25))

    # Add label outside the plot
    ax.text(-0.1, 0.95, outside_label[n], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

plt.savefig(file_out + '.png')
plt.savefig(file_out + '.pdf')
#plt.show()
