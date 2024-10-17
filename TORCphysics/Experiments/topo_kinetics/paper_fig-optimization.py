import numpy as np
import pandas as pd
from TORCphysics import topo_calibration_tools as tct
from TORCphysics import Circuit
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em
from TORCphysics import params as papa
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# Plot results from calibration process.
# And adding the binding and effect model curves

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0 #0.25
initial_time = 0
final_time = 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'topo_optimization-dt'+str(dt)
#file_out = 'calibration_small_dt1'

# For the simulation
circuit_filename = 'circuit.csv'
sites_filename = None  # 'sites_test.csv'
enzymes_filename = None  # 'enzymes_test.csv'
#environment_filename = 'environment.csv'
environment_filename = 'environment_small.csv'

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
sigma_0_topo = -0.11#0.075  # Approximately -20 supercoils according the paper
sigma_0_gyrase = 0.0  # We suppose this one.
sigma_f_gyrase = -0.11  # We also assume this one, which is the maximum at which gyrase acts.
# At this value the torque is too strong.

output_prefix = 'test0'
series = True
continuation = False

# For parallelization and calibration
n_simulations = 100

# params_file
#topo_params_file = 'calibration_dt'+str(dt)+'.csv'
#recognition_params_file = 'calibration2.csv'
recognition_params_file = 'calibration_dt'+str(dt)+'.csv'
poisson_params_file = 'calibration_'+str(dt)+'.csv'


# Models to calibrate to calibrate
# -----------------------------------
# Topoisomerase I
topoI_name = 'topoI'
topoI_type = 'environmental'
topoI_binding_model_name = 'TopoIRecognition'
topoI_effect_model_name = 'TopoILinear'
topoI_unbinding_model_name = 'PoissonUnBinding'
topoI_params = 'recognition-linear/calibration_dt'+str(dt) +'_topoI.csv'

# Gyrase
gyrase_name = 'gyrase'
gyrase_type = 'environmental'
gyrase_binding_model_name = 'GyraseRecognition'
gyrase_effect_model_name = 'GyraseLinear'
gyrase_unbinding_model_name = 'PoissonUnBinding'
gyrase_params = 'recognition-linear/calibration_dt'+str(dt) +'_gyrase.csv'

# -----------------------------------
# FIGURE Params
# -----------------------------------
width = 7
height = 3.5
lw = 4
font_size = 12
xlabel_size = 14
title_size = 16
experiment_color = 'black'
rec_model_color = 'green'
# poi_model_color = 'green'

gyrase_color = 'blue'
topoI_color = 'red'


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

    global_dict_both_sc = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                        'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                        'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                        'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma_0_topo,
                        'DNA_concentration': 0.0}

    global_dict_both_rx = {'circuit_filename': circuit_filename, 'sites_filename': sites_filename,
                        'enzymes_filename': enzymes_filename, 'environment_filename': environment_filename,
                        'output_prefix': output_prefix, 'series': series, 'continuation': continuation,
                        'frames': frames, 'dt': dt, 'n_simulations': n_simulations, 'initial_sigma': sigma_0_gyrase,
                        'DNA_concentration': 0.0}

    # Variation dictionaries
    # ------------------------------------------

    # Topoisomerase I
    name = topoI_name
    object_type = topoI_type
    binding_model_name = topoI_binding_model_name
    if 'Poisson' in binding_model_name:
        binding_oparams = {'k_on': float(params['k_on_topoI'][0])}
    else:
        binding_oparams = {'k_on': float(params['k_on_topoI'][0]), 'width': float(params['width_topoI'][0]),
                           'threshold': float(params['threshold_topoI'][0])}
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
    if 'Poisson' in binding_model_name:
        binding_oparams = {'k_on': float(params['k_on_gyrase'][0])}
    else:
        binding_oparams = {'k_on': float(params['k_on_gyrase'][0]), 'width': float(params['width_gyrase'][0]),
                           'threshold': float(params['threshold_gyrase'][0])}
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
    global_dict_list = [global_dict_topoI, global_dict_gyrase, global_dict_both_sc, global_dict_both_rx]

    # List of lists of variations
    variations_list = [[topoI_variation], [gyrase_variation],
                       [topoI_variation, gyrase_variation],
                       [topoI_variation, gyrase_variation]
                       ]

    # Arrays with global superhelical densities
    list_sigmas = [topoI_sigma, gyrase_sigma, both_sigma_sc, both_sigma_rx]

    # Finally, run objective function. run_objective_function will process our conditions
    # ------------------------------------------
    my_objective, simulation_superhelicals, simulation_superhelicals_std = tct.run_objective_function_g_std(global_dict_list=global_dict_list,
                                                                        variations_list=variations_list,
                                                                        exp_superhelicals=list_sigmas,
                                                                        n_simulations=n_simulations)
    return my_objective, simulation_superhelicals, simulation_superhelicals_std


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
# Note that this sigmaf is not at the end of the simulation, but the sigma at which there is 0 Relaxed DNA.
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
topoI_sigma = sigma

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
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
gyrase_sigma = sigma

# -----------------------------------------------------
# -----------------------------------------------------
# Build experimental curve for system with both Topo I and Gyrase (from Sc state)
# -----------------------------------------------------
# Kinetics Gyrase: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Kinetics Topoisomerase: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI

# Integrate MM kinetics
# ------------------------------------------
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_both_T_G(vmax_topoI=v_max_topoI, vmax_gyrase=v_max_gyrase,
                                                         KM_topoI=K_M_topoI, KM_gyrase=K_M_gyrase,
                                                         Supercoiled_0=DNA_concentration, Relaxed_0=0.0,
                                                         frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
both_sigma_sc = sigma

# -----------------------------------------------------
# Build experimental curve for system with both Topo I and Gyrase (from Rx state)
# -----------------------------------------------------
# Integrate MM kinetics
# ------------------------------------------
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_both_T_G(vmax_topoI=v_max_topoI, vmax_gyrase=v_max_gyrase,
                                                         KM_topoI=K_M_topoI, KM_gyrase=K_M_gyrase,
                                                         Supercoiled_0=0.0, Relaxed_0=DNA_concentration,
                                                         frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
both_sigma_rx = sigma

# ----------------------------------------------------------------------------------------------------------------------
# Theoretical curves
# ----------------------------------------------------------------------------------------------------------------------

# Recognition binding
#params_file = 'recognition-linear/calibration.csv'
#params_file = 'recognition-linear/calibration_small_dt1.csv'
params_file = 'recognition-linear/'+recognition_params_file
params_dict = pd.read_csv(params_file).to_dict()
rec_objective, rec_sim_superhelicals, rec_sim_superhelicals_std = objective_function(params=params_dict)

# ----------------------------------------------------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------------------------------------------------

# Create figure
# fig, axs = plt.subplots(2, 2, figsize=(width*2, 2 * height), tight_layout=True)
fig, axs = plt.subplots(3, 2, figsize=(width*2, 3 * height), tight_layout=True)

# Prepare arrays to iterate
exp_sigmas = [topoI_sigma, gyrase_sigma, both_sigma_sc, both_sigma_rx]
titles = ['Topoisomerase I Acting on Supercoiled (SC) DNA', 'Gyrase Acting on Relaxed (RX) DNA',
          'Topoisomerase I & Gyrase Acting on SC DNA', 'Topoisomerase I & Gyrase Acting on RX DNA',
          'Topoisomerase Recognition Binding Model', 'Topoisomerase Linear Effect Model']

outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']
err_x = [0.75, 0.75, 0.75, 0.75]
err_y = [0.4, 0.4 ,0.35, 0.35]

# Plot calibration
# ------------------------------------------
axes = [axs[0,0], axs[0,1], axs[1,0], axs[1,1]]
ylim2 = 0.01
ylim1 = -0.12
for n in range(4):
    ax = axes[n]
    sum_err_rec = np.sum(np.square(rec_sim_superhelicals[n] - exp_sigmas[n]))/len(rec_sim_superhelicals[n])

    ax.plot(time, exp_sigmas[n], lw=lw, color=experiment_color, label='Kinetic')
    ax.plot(time, rec_sim_superhelicals[n], lw=lw, color=rec_model_color, label='TORCPhysics')
    y = rec_sim_superhelicals[n]
    ys =rec_sim_superhelicals_std[n]
    ax.fill_between(time, y-ys, y+ys, lw=lw, color=rec_model_color, alpha=0.2)

    ax.set_xlabel('Time (s)', fontsize=xlabel_size)
    ax.set_ylabel('Superhelical density', fontsize=xlabel_size)
    ax.grid(True)
    ax.set_title(titles[n], fontsize=title_size)
    ax.set_ylim(ylim1, ylim2)

    # Add errors
    formatted_sum_err = "{:.3g}".format(sum_err_rec)  # Formatting to three significant figures
    er_label = r'$\epsilon={}$'.format(formatted_sum_err)
    props = dict(boxstyle='round', facecolor='silver', alpha=0.4)
    ax.text(err_x[n], err_y[n], er_label, transform=ax.transAxes, fontsize=font_size,
            verticalalignment='top', bbox=props)


    # Add label outside the plot
    ax.text(-0.07, 1.0, outside_label[n], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

axs[0,0].legend(loc='best', fontsize=font_size)

# Plot Models
# --------------------------------------------------------------------------
sigma = np.arange(-.15, .15, 0.01)

# Initialize circuit
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt)

# Responses
# -------------------------
ax = axs[2,0]
n = 4
gyrase_bind_rate = gyrase_concentration * bm.GyraseRecognition(filename=gyrase_params).rate_modulation(superhelical=sigma)
topoI_bind_rate = topoI_concentration * bm.TopoIRecognition(filename=topoI_params).rate_modulation(superhelical=sigma)

ax.plot(sigma, gyrase_bind_rate, color=gyrase_color, lw=lw, label='Gyrase')
ax.plot(sigma, topoI_bind_rate, color=topoI_color, lw=lw, label='Topoisomerase I')

ax.set_title(titles[n], fontsize=title_size)
ax.set_xlabel('Superhelical Density', fontsize=xlabel_size)
ax.set_ylabel(r'Binding Rate ($s^{-1}$)', fontsize=xlabel_size)
# Add label outside the plot
ax.text(-0.07, 1.0, outside_label[n], transform=ax.transAxes,
        fontsize=font_size * 1.5, fontweight='bold', va='center', ha='center')
ax.legend(loc='best', fontsize=font_size)
ax.grid(True)

# Effects
# -------------------------
ax = axs[2,1]
n=5
gyrase_ef = em.GyraseLinear(filename=gyrase_params).twist_added(superhelical=sigma, dt=1.0)/papa.w0
topoI_ef = em.TopoILinear(filename=topoI_params).twist_added(superhelical=sigma, dt=1.0)/papa.w0
ax.plot(sigma, gyrase_ef, color=gyrase_color, lw=lw)
ax.plot(sigma, topoI_ef, color=topoI_color, lw=lw)

ax.set_title(titles[n], fontsize=title_size)
ax.set_xlabel('Superhelical Density', fontsize=xlabel_size)
ax.set_ylabel(r'Twist Induced (bp)', fontsize=xlabel_size)
# Add label outside the plot
ax.text(-0.07, 1.0, outside_label[n], transform=ax.transAxes,
        fontsize=font_size * 1.5, fontweight='bold', va='center', ha='center')
#ax.legend(loc='best', fontsize=font_size)
ax.grid(True)

plt.savefig(file_out + '.png')
plt.savefig(file_out + '.pdf')
plt.show()
