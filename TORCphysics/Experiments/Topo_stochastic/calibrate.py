from TORCphysics import Circuit
import matplotlib.pyplot as plt
from TORCphysics import visualization as vs
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import sys
import random
import numpy as np
from scipy.optimize import curve_fit

nsim = 3


# Here, we want to focuse on k_cat
def stochastic_topo_kcat(time, topo_kcat, gyra_kcat):
    # change params
    nframes = len(time) - 1
    supercoiling = np.zeros((nframes + 1, nsim))

    for n in range(nsim):
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, nframes, series, continuation, dt, tm, mm)
        my_circuit.environmental_list[0].k_on = topo_kon0
        my_circuit.environmental_list[0].k_cat = topo_kcat
        my_circuit.environmental_list[1].k_on = gyra_kon0
        my_circuit.environmental_list[1].k_cat = gyra_kcat

        supercoiling[0, n] = my_circuit.superhelical

        # run simulation
        for frame in range(1, nframes + 1):
            # BINDING
            # --------------------------------------------------------------
            new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                               my_circuit.rng)
            my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

            # EFFECT
            # --------------------------------------------------------------
            effects_list = em.effect_model(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                           my_circuit.topoisomerase_model, my_circuit.mechanical_model)
            my_circuit.apply_effects(effects_list)

            # UNBINDING
            # --------------------------------------------------------------
            drop_list_index, drop_list_enzyme = bm.unbinding_model(my_circuit.enzyme_list, my_circuit.dt,
                                                                   my_circuit.rng)
            my_circuit.drop_enzymes(drop_list_index)
            my_circuit.add_to_environment(drop_list_enzyme)
            # UPDATE GLOBALS
            # --------------------------------------------------------------
            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            supercoiling[frame, n] = my_circuit.superhelical
    return np.mean(supercoiling, axis=1)


# This function runs a simulation for calibrating the topo activity
def stochastic_topo(time, topo_kon, topo_kcat, gyra_kon, gyra_kcat):
    # change params
    nframes = len(time) - 1
    supercoiling = np.zeros((nframes + 1, nsim))

    for n in range(nsim):
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, nframes, series, continuation, dt, tm, mm)
        my_circuit.environmental_list[0].k_on = topo_kon
        my_circuit.environmental_list[0].k_cat = topo_kcat
        my_circuit.environmental_list[1].k_on = gyra_kon
        my_circuit.environmental_list[1].k_cat = gyra_kcat

        supercoiling[0, n] = my_circuit.superhelical

        # run simulation
        for frame in range(1, nframes + 1):
            # BINDING
            # --------------------------------------------------------------
            new_enzyme_list = bm.binding_model(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                               my_circuit.rng)
            my_circuit.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

            # EFFECT
            # --------------------------------------------------------------
            effects_list = em.effect_model(my_circuit.enzyme_list, my_circuit.environmental_list, my_circuit.dt,
                                           my_circuit.topoisomerase_model, my_circuit.mechanical_model)
            my_circuit.apply_effects(effects_list)

            # UNBINDING
            # --------------------------------------------------------------
            drop_list_index, drop_list_enzyme = bm.unbinding_model(my_circuit.enzyme_list, my_circuit.dt,
                                                                   my_circuit.rng)
            my_circuit.drop_enzymes(drop_list_index)
            my_circuit.add_to_environment(drop_list_enzyme)
            # UPDATE GLOBALS
            # --------------------------------------------------------------
            my_circuit.update_global_twist()
            my_circuit.update_global_superhelical()
            supercoiling[frame, n] = my_circuit.superhelical
    return np.mean(supercoiling, axis=1)


# TODO:
#  1.- Make the code work with stochastic topo binding. DONE
#  2.- Run experiment to test that it works. DONE
#  3.- Calibrate it so it works with Sam Meyers model, for now...
#  3.1.- I need to plot the global supercoiling response of the continuum case
#  3.2.- Plot the response of my initial case
#  3.3.- Run multiple experiments using the curve fit to find the optimal parameters of my stochastic topo binding

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# I will be plotting as the experiment advances.

# Initialise figure and circuit initial conditions
# ----------------------------------------------------------------------------------------------------------------------

# Figure
width = 8
height = 3
figure_output = 'calibration.png'

fig, axs = plt.subplots(2, figsize=(width, 2.5 * height), tight_layout=True)

# Circuit conditions
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_filename = 'environment_continuum.csv'
output_prefix = 'continuum'
frames = 750
series = True
continuation = False
mm = 'uniform'
dt = 1.0

# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------
tm = 'continuum'
continuum_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                            output_prefix, frames, series, continuation, dt, tm, mm)
continuum_circuit.run()

# Get global supercoiling responses from the continuum case
mask = continuum_circuit.sites_df['type'] == 'circuit'  # This one contains global superhelical density
sigma_continuum = continuum_circuit.sites_df[mask]['superhelical'].to_numpy()

# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------
tm = 'stochastic'
output_prefix = 'stochastic'
environment_filename = 'environment.csv'
# Load simulation
test_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                       output_prefix, frames, series, continuation, dt, tm, mm)

# Let's check if the function works
time = np.arange(0, (frames + 1) * dt, dt)
topo_kon0 = .005  # 0.0075
topo_kcat0 = 20
gyra_kon0 = .005  # 0.0075
gyra_kcat0 = -20
test_circuit.environmental_list[0].k_on = topo_kon0
test_circuit.environmental_list[0].k_cat = topo_kcat0
test_circuit.environmental_list[1].k_on = gyra_kon0
test_circuit.environmental_list[1].k_cat = gyra_kcat0

p0 = [topo_kon0, topo_kcat0, gyra_kon0, gyra_kcat0]
pcat0 = [1, -1]
sigma_stochastic0 = stochastic_topo(time, *p0)
ax = axs[0]
ax.plot(time, sigma_continuum, color='black')
ax.plot(time, sigma_stochastic0, color='blue')
# sigma_stochastic2 = stochastic_topo(time, topo_kon=0.003, topo_kcat=3, gyra_kon=0.003, gyra_kcat=-3)
# ax.plot(time, sigma_stochastic2, color='green')
# popt, pcov = curve_fit(stochastic_topo, time, sigma_continuum, #  method='dogbox',
#                       p0=p0, bounds=((0.0001, 1, 0.0001, -50), (0.1, 50, 0.1, -1)))
# sigma_stochastic = stochastic_topo(time, *popt)
popt, pcov = curve_fit(stochastic_topo_kcat, time, sigma_continuum,  # method='dogbox',
                       p0=pcat0)  #), bounds=((1, 0.0001, -50), (0.1, 50, 0.1, -1)))
sigma_stochastic = stochastic_topo_kcat(time, *popt)

ax.plot(time, sigma_stochastic, color='red')
print(popt)
ax.grid(True)
# n_simulations = 1
# for ns in range(n_simulations):

#    # Load simulation
#    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
#                         output_prefix, frames, series, continuation, dt, tm, mm)

#    my_circuit.name = my_circuit.name + '_' + str(ns)
#    my_circuit.log.name = my_circuit.name
#  my_circuit.print_general_information()
#    my_circuit.run()


# Plot site response curves
# ---------------------------------------------------------
ax = axs[1]
vs.plot_site_response_curves(test_circuit, ax)  # , ignore=to_ignore)
# Let's ignore the

# Plot global supercoiling responses
# ---------------------------------------------------------
# ax = axs[2]
# vs.plot_supercoiling_profiles(my_circuit, my_circuit.sites_df, ax, only_global=True)

plt.savefig('calibration.png')
