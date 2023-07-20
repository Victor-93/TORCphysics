import matplotlib.pyplot as plt
import numpy as np
from TORCphysics import Circuit

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This is just a test to reproduce the global supercoiling response curves from the paper:
# Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching
# Here, I want to test how can I use,

# TODO:
#  1.- Do it for gyrase response.
#  2.- See how you set the initial supercoiling and final (probably it does not tend to 0).
#  3.- See how you can add the kinetics with ATP to your model.
#  4.- And check weather you need to update the shape of your recognition curves (probably not, or hopefully).
#  5.- And finally, check if from the experimental K_M, you can propose your k_on and k_off.

# TODO nuevo:
# 1- Vamos a juntar con las curvas de Meyer, con la misma concentracion y tiempo de simulacion, a ver que pasa
# 2- Luego vamos tenemos que definir la sigma inicial y final.
# 3 - Una vez haciendo eso, haces lo de arriba en ingles

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1
initial_time = 0
final_time = 600
initial_sigma = -.04
final_sigma = 0.0
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

# Meyer topos curve parameters:
topo_w = 0.012  # width
topo_t = -0.04  # thresholds
topo_k = 0.001  # k_cat
gyra_w = 0.025  # width
gyra_t = 0.01  # threshold
gyra_k = 0.001  # k_cat

# Figure
width = 6
height = 3
lw = 3
exp_color = 'blue'
sam_color = 'red'
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
initial_product = 0.0
initial_substrate = .7
enzyme_concentration = 17
K_M = 1.5
k_cat = .0023  # 0.003
v_max = k_cat * enzyme_concentration
# Substrates and products
# ----------------------------------
ax = axs[0, 0]
substrate, product = integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                  frames=frames, dt=dt)
ax.plot(time, product, lw=lw, color=exp_color)

# Sigma deduction
# ----------------------------------
ax = axs[1, 0]
superhelical = rescale_product_to_sigma(product, initial_sigma, final_sigma)
ax.plot(time, superhelical, lw=lw, color=exp_color)

# Let's run the sam Meyer simulation
# ----------------------------------
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt, tm, mm)
my_circuit.environmental_list[0].concentration = enzyme_concentration/1000  # topoI
my_circuit.environmental_list[0].k_cat = 0.02 #topo_k  # topoI
my_circuit.environmental_list[1].concentration = 0.0  # gyrase
my_circuit.superhelical = initial_sigma
my_circuit.run()

mask = my_circuit.sites_df['type'] == 'circuit'  # This one contains global superhelical density
superhelical_continuum = my_circuit.sites_df[mask]['superhelical'].to_numpy()

ax.plot(time, superhelical_continuum[0:frames], lw=lw, color=sam_color)

# Gyrase
# ==================================================================================================================
# Kinetics: RDNA + Gyrase + ATP -> RDNA-Gyrase-ATP -> SDNA + Gyrase + ADP
# Product = Less fluorescent or Supercoiled DNA
# Substrate = Concentration of Relaxed DNA or fluorescent
initial_product = 4.0
initial_substrate = .75
enzyme_concentration = 44.6
K_M = 2.7
k_cat = .0011
v_max = k_cat * enzyme_concentration
# Substrates and products
# ----------------------------------
ax = axs[0, 1]
substrate, product = integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                  frames=frames, dt=dt)
ax.plot(time, substrate, lw=lw)
# Sigma deduction
# ----------------------------------
ax = axs[1, 1]
superhelical = rescale_product_to_sigma(substrate, initial_sigma, final_sigma)
ax.plot(time, superhelical, lw=lw)

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
