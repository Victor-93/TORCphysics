import matplotlib.pyplot as plt
import numpy as np
from TORCphysics import Circuit
from hyperopt import tpe, hp, fmin
import multiprocessing

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This is just a test to reproduce the global supercoiling response curves from the paper:
# Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching

# We now have estimations of the initial and final superhelical densities.
# For now, let's not consider ATP.
# Let's run examples of the actual calibration. Remember that in your simulation, you need to consider the
# density of the plasmid concentration as well, so you have many concentrations to test.

# TODO:
#  1.- Run parallelization test, varying plasmid concentrations.
#  2.- MY initial thought was that according K_M & k_cat, I could determine k_on and k_off, but probably
#  I'll need to determine these parameters. I think I can consider a rapid disassociation, and slow binding.
#  3.- So in general, the parameters I need to optimize are k_on, k_off, k_cat, width and threshold.
# TODO: Remember that you need to add the substrate concentration


# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1
initial_time = 0
final_time = 600
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# For the simulation
circuit_filename = 'circuit_test.csv'
sites_filename = 'sites_test.csv'
enzymes_filename = 'enzymes_test.csv'
environment_filename = 'environment_test.csv'
tm = 'continuum'
output_prefix = 'test0'
series = True
continuation = False
mm = 'uniform'
dt_sim = 0.5

# For parallelization and calibration
n_simulations = 5
tests = 3 #00  # number of tests for parametrization

my_vars = ['k_cat', 'k_on', 'k_off', 'width', 'threshold']

# Meyer topos curve parameters:
topo_w = 0.012  # width
topo_t = -0.04  # thresholds
topo_k = 0.001  # k_cat
gyra_w = 0.025  # width
gyra_t = 0.01  # threshold
gyra_k = 0.001  # k_cat

# Gyrase range
coutput = 'gyrase_calibration.test'
k_cat_min = -20.0  # Ranges to vary k_cat
k_cat_max = -5.0
width_min = 0.001
width_max = 5.
threshold_min = 0.001
threshold_max = 5.

k_off = 0.5
k_on = 0.005
concentration = 0.25


# Figure
width = 6
height = 3
lw = 3
colors = ['blue', 'red', 'green', 'black', 'grey', 'yellow', 'cyan', 'purple', 'brown', 'pink']
fig, axs = plt.subplots(2, 2, figsize=(2 * width, 2 * height), tight_layout=True)


# ----------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------
def Michael_Menten_equation(vmax, KM, S):
    return vmax * S / (KM + S)


def integrate_MM(vmax, KM, substrate0, product0, frames, dt):
    substrate_array = np.zeros(frames)
    product_array = np.zeros(frames)
    substrate_array[0] = substrate0
    product_array[0] = product0

    substrate = substrate0
    product = product0

    for k in range(1, frames):
        v = Michael_Menten_equation(vmax=vmax, KM=KM, S=substrate)
        product = product + v * dt
        substrate = substrate - v * dt
        substrate_array[k] = substrate
        product_array[k] = product
    return substrate_array, product_array


def rescale_product_to_sigma(product, sigma_min, sigma_max):
    #    current_min
    #    frames = len(product)
    #    sigma = np.zeros(frames)
    #    sigma = product-np.min(product)
    #    sigma = sigma/np.max(sigma) # Normalized
    #    d = sigma_f - sigma_i
    #    sigma = (sigma - sigma_i)/1#abs(sigma_f-sigma_i)

    current_min = np.min(product)
    current_max = np.max(product)
    range_ = current_max - current_min
    desired_range = sigma_max - sigma_min
    sigma = ((product - current_min) / range_) * desired_range + sigma_min
    return sigma


# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------

# TOPOISOMERASE I
# ==================================================================================================================
# Kinetics: SDNA + TopoI -> SDNA-TopoI -> RDNA + TopoI
# Product = Fluorescent or Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs
initial_sigma = -.047
final_sigma = 0.0
initial_product = 0.0
initial_substrate = .7
initial_substrates = [0.35, 0.7, 1.1, 1.8, 2.1]
enzyme_concentration = 17
K_M = 1.5
k_cat = .0023  # 0.003
v_max = k_cat * enzyme_concentration
for count, initial_substrate in enumerate(initial_substrates):
    # Substrates and products
    # ----------------------------------
    ax = axs[0, 0]
    substrate, product = integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                      frames=frames, dt=dt)
    ax.plot(time, product, lw=lw, color=colors[count])

    # Sigma deduction
    # ----------------------------------
    ax = axs[1, 0]
    superhelical = rescale_product_to_sigma(product, initial_sigma, final_sigma)
    ax.plot(time, superhelical, lw=lw, color=colors[count])

# Gyrase
# ==================================================================================================================
# Kinetics: RDNA + Gyrase + ATP -> RDNA-Gyrase-ATP -> SDNA + Gyrase + ADP
# Product = Less fluorescent or Supercoiled DNA
# Substrate = Concentration of Relaxed DNA or fluorescent
initial_sigma = -0.02
final_sigma = 0.0
initial_product = 4.0
initial_substrate = .75
initial_substrates = [0.75, 1.50, 3.6, 5.4, 7.2]
enzyme_concentration = 44.6
K_M = 2.7
k_cat = .0011
v_max = k_cat * enzyme_concentration
for count, initial_substrate in enumerate(initial_substrates):
    # Substrates and products
    # ----------------------------------
    ax = axs[0, 1]
    substrate, product = integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                      frames=frames, dt=dt)
    ax.plot(time, substrate, lw=lw, color=colors[count])

    # Sigma deduction
    # ----------------------------------
    ax = axs[1, 1]
    superhelical = rescale_product_to_sigma(substrate, initial_sigma, final_sigma)
    ax.plot(time, superhelical, lw=lw, color=colors[count])

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

plt.show()
