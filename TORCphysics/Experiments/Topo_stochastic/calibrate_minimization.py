from TORCphysics import Circuit
import matplotlib.pyplot as plt
from TORCphysics import visualization as vs
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import sys
import random
import numpy as np
from scipy.optimize import curve_fit, minimize

# Now we will calibrate by minimizing an objective function

nsim = 20
# Circuit conditions
circuit_filename = 'circuit.csv'
sites_filename = 'sites.csv'
enzymes_filename = 'enzymes.csv'
environment_continuum_filename = 'environment_continuum.csv'
environment_filename = 'environment.csv'
tm = 'stochastic'
output_prefix = 'continuum'
nframes = 500
series = True
continuation = False
mm = 'uniform'
dt = 1.0


# Here, we want to focuse on k_cat

# This function runs a simulation for calibrating the topo activity
def objective_fun_1(x):
    # Run continuum case
    # ----------------------------------------------------------------------------------------------------------------------
    tm = 'continuum'
    continuum_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_continuum_filename,
                                output_prefix, nframes, series, continuation, dt, 'continuum', mm)
    continuum_circuit.run()
    # Get global supercoiling responses from the continuum case
    mask = continuum_circuit.sites_df['type'] == 'circuit'  # This one contains global superhelical density
    sigma_continuum = continuum_circuit.sites_df[mask]['superhelical'].to_numpy()

    # change params
    supercoiling = np.zeros((nframes + 1, nsim))
    for n in range(nsim):
        supercoiling[:, n] = run_stochastic_sim(x, circuit_filename, sites_filename, enzymes_filename,
                                                environment_filename, output_prefix, nframes,
                                                series, continuation, dt, mm)
    objective = np.sum(np.square(np.mean(supercoiling, axis=1) - sigma_continuum))
    return objective


def run_stochastic_sim(x, circuit_filename, sites_filename, enzymes_filename, environment_filename,
                       output_prefix, nframes, series, continuation, dt, mm):
    my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                         output_prefix, nframes, series, continuation, dt, 'stochastic', mm)
    my_circuit.environmental_list[0].k_on = x[0]
    my_circuit.environmental_list[0].k_cat = x[1]
    my_circuit.environmental_list[1].k_on = x[2]
    my_circuit.environmental_list[1].k_cat = x[3]

    my_supercoiling = np.zeros(nframes + 1)
    my_supercoiling[0] = my_circuit.superhelical

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
        my_supercoiling[frame] = my_circuit.superhelical
    return my_supercoiling

def objective_function_kcat(topo_kcat, gyra_kcat, dcat):
    tm = 'continuum'
    continuum_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_continuum_filename,
                                output_prefix, nframes, series, continuation, dt, 'continuum', mm)
    continuum_circuit.run()
#    # Get global supercoiling responses from the continuum case
    mask = continuum_circuit.sites_df['type'] == 'circuit'  # This one contains global superhelical density
    sigma_continuum = continuum_circuit.sites_df[mask]['superhelical'].to_numpy()

    topo_kcat_ranges = np.arange(topo_kcat[0], topo_kcat[1]+dcat, dcat)
    gyra_kcat_ranges = np.arange(gyra_kcat[0], gyra_kcat[1]+dcat, dcat)

    kcat_list = []
    objective = []
    #objective = np.ones((len(topo_kcat_ranges), len(gyra_kcat_ranges)))*10
    for i, topocat in enumerate(topo_kcat_ranges):
        for j, gyracat in enumerate(gyra_kcat_ranges):
            kcat_list.append( [topocat, gyracat] )
            x = [0.005, topocat, 0.005, gyracat]
            supercoiling = np.zeros((nframes + 1, nsim))
            for n in range(nsim):
                supercoiling[:, n] = run_stochastic_sim(x, circuit_filename, sites_filename, enzymes_filename,
                                                        environment_filename, output_prefix, nframes,
                                                        series, continuation, dt, mm)
            objective.append( np.sum(np.square(np.mean(supercoiling, axis=1) - sigma_continuum)) )
    return objective, kcat_list



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
figure_output = 'calibration_min.png'

fig, axs = plt.subplots(2, figsize=(width, 2.5 * height), tight_layout=True)

time = np.arange(0, (nframes + 1) * dt, dt)

# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------
tm = 'continuum'
continuum_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_continuum_filename,
                            output_prefix, nframes, series, continuation, dt, tm, mm)
continuum_circuit.run()

# Get global supercoiling responses from the continuum case
mask = continuum_circuit.sites_df['type'] == 'circuit'  # This one contains global superhelical density
sigma_continuum = continuum_circuit.sites_df[mask]['superhelical'].to_numpy()

ax = axs[0]
ax.plot(time, sigma_continuum, color='black')

# Run continuum case
# ----------------------------------------------------------------------------------------------------------------------

# Other case
# ----------------------------------------------------------------------------------------------------------------------
tcats = [2.5, 30]
gcats = [-20,-2.5]
objective, kcat_list = objective_function_kcat(tcats, gcats, dcat=2.5)
objective = np.array(objective)
print(objective)
print(np.min(objective))
a = np.argmin(objective)
print(kcat_list[a])
x0 = [0.005, kcat_list[a][0], 0.005, kcat_list[a][1]]
sigma_stochastic_0 = np.zeros((len(sigma_continuum), nsim))
sigma_stochastic = np.zeros((len(sigma_continuum), nsim))
for n in range(nsim):
#    sigma_stochastic[:, n] = run_stochastic_sim(res.x, circuit_filename, sites_filename, enzymes_filename,
#                                                environment_filename, output_prefix, nframes, series, continuation, dt,
#                                                mm)
    sigma_stochastic_0[:, n] = run_stochastic_sim(x0, circuit_filename, sites_filename, enzymes_filename,
                                                  environment_filename,
                                                  output_prefix, nframes, series, continuation, dt, mm)

ax.plot(time, np.mean(sigma_stochastic_0, axis=1), color='blue')
ax.grid(True)

plt.savefig(figure_output)

sys.exit()

# Let's check if the function works
topo_kon0 = .005  # 0.0075
topo_kcat0 = 15
gyra_kon0 = .005  # 0.0075
gyra_kcat0 = -0
x0 = [topo_kon0, topo_kcat0, gyra_kon0, gyra_kcat0]



#res = minimize(objective_fun_1, x0, method='BFGS')  # ), bounds=((1, 0.0001, -50), (0.1, 50, 0.1, -1)))
res = minimize(objective_fun_1, x0, method='Nelder-Mead', tol=1e-6)
print(res)
print(res.x)


sigma_stochastic_0 = np.zeros((len(sigma_continuum), nsim))
sigma_stochastic = np.zeros((len(sigma_continuum), nsim))
for n in range(nsim):
    sigma_stochastic[:, n] = run_stochastic_sim(res.x, circuit_filename, sites_filename, enzymes_filename,
                                                environment_filename, output_prefix, nframes, series, continuation, dt,
                                                mm)
    sigma_stochastic_0[:, n] = run_stochastic_sim(x0, circuit_filename, sites_filename, enzymes_filename,
                                                  environment_filename,
                                                  output_prefix, nframes, series, continuation, dt, mm)

objective = np.sum(np.square(np.mean(sigma_stochastic, axis=1) - sigma_continuum))
print(objective)
ax.plot(time, np.mean(sigma_stochastic_0, axis=1), color='blue')
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
# vs.plot_site_response_curves(test_circuit, ax)  # , ignore=to_ignore)
# Let's ignore the

# Plot global supercoiling responses
# ---------------------------------------------------------
# ax = axs[2]
# vs.plot_supercoiling_profiles(my_circuit, my_circuit.sites_df, ax, only_global=True)

plt.savefig(figure_output)
