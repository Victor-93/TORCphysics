import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This is just a test to reproduce the global supercoiling response curves from the paper:
# Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
initial_substrate = .35
initial_substrates = [2.1]  # [2.1, 1.8, 1.1, .7, .35]
dt = 1
initial_time = 0
final_time = 400
initial_sigma = -.04
final_sigma = 0.0
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)

# Figure
width = 6
height = 3
fig, axs = plt.subplots(2, 2, figsize=(1.25 * width, 2 * height), tight_layout=True)


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

for initial_substrate in initial_substrates:
    # TOPOISOMERASE I
    # ==================================================================================================================
    initial_product = 0.0
    enzyme_concentration = 17
    K_M = 1.5
    k_cat = .0023  # 0.003
    v_max = k_cat * enzyme_concentration
    # Substrates and products
    # ----------------------------------
    ax = axs[0, 0]
    substrate, product = integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                      frames=frames, dt=dt)
    ax.plot(time, product)
    ax.plot(time, substrate)

    # Sigma deduction
    # ----------------------------------
    ax = axs[1, 0]
    superhelical = rescale_product_to_sigma(product, initial_sigma, final_sigma)
    ax.plot(time, superhelical)

    # Gyrase
    # ==================================================================================================================
    initial_product = 4.0
    enzyme_concentration = 44.6
    K_M = 2.7
    k_cat = .0011
    v_max = -k_cat * enzyme_concentration
    # Substrates and products
    # ----------------------------------
    ax = axs[0, 1]
    substrate, product = integrate_MM(vmax=v_max, KM=K_M, substrate0=initial_substrate, product0=initial_product,
                                      frames=frames, dt=dt)
    ax.plot(time, product)
    # Sigma deduction
    # ----------------------------------
    ax = axs[1, 1]
    superhelical = rescale_product_to_sigma(product, initial_sigma, final_sigma)
    ax.plot(time, superhelical)

for ax in [axs[0,0], axs[0,1]]:
    ax.set_xlabel('time (s)')
    ax.set_ylabel('Relaxed DNA (nM)')
    ax.grid(True)

for ax in [axs[1,0], axs[1,1]]:
    ax.set_xlabel('time (s)')
    ax.set_ylabel('Superhelical response')
    ax.grid(True)

axs[0, 0].set_title('Topoisomerase I')
axs[0, 1].set_title('Gyrase')

plt.show()
