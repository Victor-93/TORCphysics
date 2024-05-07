import numpy as np
import matplotlib.pyplot as plt
from TORCphysics import topo_calibration_tools as tct

# ----------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ----------------------------------------------------------------------------------------------------------------------
# This script is to make sure that our reproduced kinetic curves are in agreement with literature.
# We want to plot the kinetic curves for Topo I, Gyrase and both enzymes acting together.
# For each case, we will produce three quantities as a function of time:
# 1.- [P] = [R] Product or Relaxed DNA concentration.
# 2.- [S] Substrate or Supercoiled DNA.
# 3.- sigma - superhelical density

# We want to plot the results of the calibration

# ----------------------------------------------------------------------------------------------------------------------
# Initial conditions
# ----------------------------------------------------------------------------------------------------------------------
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 0.25
initial_time = 0
final_time = 400
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'kinetic_curves'

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
sigma_0_topo = -0.11 #-0.076  # Approximately -20 supercoils according the paper
sigma_0_gyrase = 0.0  # We suppose this one.
sigma_f_gyrase = -0.11  # We also assume this one, which is the maximum at which gyrase acts.
# At this value the torque is too strong.

# -----------------------------------
# FIGURE
# -----------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16
topo_color = 'red'
gyrase_color = 'blue'
both_color = 'green'
fig, axs = plt.subplots(3,  figsize=(width, 3 * height), tight_layout=True)

# ----------------------------------------------------------------------------------------------------------------------
# Process
# ----------------------------------------------------------------------------------------------------------------------

# ==================================================================================================================
# Build experimental curve for TOPO I
# ==================================================================================================================
# Kinetics: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI
# Product = Relaxed DNA
# Substrate = Concentration of Supercoiled DNAs; which initially is the same as the DNA concentration
mycolor = topo_color
mylabel = 'Topo I'

# Integrate MM kinetics
# ------------------------------------------
# Initially, there's no relaxed DNA, and all the supercoiled DNA concentration corresponds to the plasmid conc.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_topoI(vmax=v_max_topoI, KM=K_M_topoI,
                                                      Supercoiled_0=DNA_concentration, Relaxed_0=0.0,
                                                      frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
#sigma = tct.topoI_to_sigma(Relaxed=relaxed_DNA, DNA_concentration=DNA_concentration, sigma0=sigma_0_topo)
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)

# Plot data
# ------------------------------------------
axs[0].plot(time, relaxed_DNA, lw=lw, color=mycolor, label=mylabel)
axs[1].plot(time, supercoiled_DNA, lw=lw, color=mycolor, label=mylabel)
axs[2].plot(time, sigma, lw=lw, color=mycolor, label=mylabel)

# ==================================================================================================================
# Build experimental curve for Gyrase
# ==================================================================================================================
# Kinetics: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Product = Supercoiled DNA
# Substrate = Relaxed DNA; which initially is the same as the DNA concentration
mycolor = gyrase_color
mylabel = 'Gyrase'

# Integrate MM kinetics
# ------------------------------------------
# Initially, there's no supercoiled DNA, and all of the relaxed DNA concentration corresponds
# to the plasmid concentration.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_gyrase(vmax=v_max_gyrase, KM=K_M_gyrase,
                                                       Supercoiled_0=0.0, Relaxed_0=DNA_concentration,
                                                       frames=frames, dt=dt)
# Translate to superhelical density
# ------------------------------------------
#sigma = tct.gyrase_to_sigma(Relaxed=relaxed_DNA, DNA_concentration=DNA_concentration,
#                            sigma0=sigma_0_gyrase, sigmaf=sigma_f_gyrase)
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)

# Plot data
# ------------------------------------------
axs[0].plot(time, relaxed_DNA, lw=lw, color=mycolor, label=mylabel)
axs[1].plot(time, supercoiled_DNA, lw=lw, color=mycolor, label=mylabel)
axs[2].plot(time, sigma, lw=lw, color=mycolor, label=mylabel)

# ==================================================================================================================
# Build experimental curve for Gyrase and topo I acting on the DNA
# ==================================================================================================================
# Kinetics Gyrase: Relaxed_DNA + Gyrase -> Relaxed-Gyrase -> Supercoiled_DNA + Gyrase
# Kinetics Topoisomerase: Supercoiled_DNA + TopoI -> Supercoiled_DNA-TopoI -> Relaxed_DNA + TopoI

mycolor = both_color
mylabel = 'Both'

# Integrate MM kinetics
# ------------------------------------------

# Initially, there's no supercoiled DNA, and all of the relaxed DNA concentration corresponds
# to the plasmid concentration.
supercoiled_DNA, relaxed_DNA = tct.integrate_MM_both_T_G(vmax_topoI=v_max_topoI, vmax_gyrase=v_max_gyrase,
                                                         KM_topoI=K_M_topoI, KM_gyrase=K_M_gyrase,
                                                         Supercoiled_0=DNA_concentration, Relaxed_0=0.0,
                                                         frames=frames, dt=dt)

ratio = relaxed_DNA[-1]/DNA_concentration
sigmaf = sigma_0_topo*ratio
sigmaf = sigma_f_gyrase - sigma_f_gyrase * relaxed_DNA[-1]/DNA_concentration
print(sigmaf)

# Translate to superhelical density
# ------------------------------------------
sigma = tct.sigma_to_relaxed(Relaxed=relaxed_DNA,
                             DNA_concentration=DNA_concentration,
                             sigmaf=sigma_f_gyrase)
#sigma = tct.both_T_G_to_sigma(Relaxed=relaxed_DNA, Relaxed_final=relaxed_DNA[-1],
#                              sigma0=sigma_0_topo, sigmaf=sigmaf)

# Plot data
# ------------------------------------------
axs[0].plot(time, relaxed_DNA, lw=lw, color=mycolor, label=mylabel)
axs[1].plot(time, supercoiled_DNA, lw=lw, color=mycolor, label=mylabel)
axs[2].plot(time, sigma, lw=lw, color=mycolor, label=mylabel)


# Sort labels
# ==================================================================================================================
outside_label = ['a)', 'b)', 'c)']
for i in range(3):

    axs[i].grid(True)
    axs[i].set_xlabel('Time (s)', fontsize=xlabel_size)

    # Add label outside the plot
    axs[i].text(-0.12, 0.99, outside_label[i], transform=axs[i].transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')


axs[0].set_ylabel('Relaxed DNA (nM)', fontsize=xlabel_size)
axs[1].set_ylabel('Supercoiled DNA (nM)', fontsize=xlabel_size)
axs[2].set_ylabel(r'Superhelical density', fontsize=xlabel_size)

axs[0].set_ylim(-.05,.85)
axs[1].set_ylim(-.05,.85)
axs[2].set_ylim(-.12,.02)
axs[2].legend(loc='best', fontsize=font_size)
    #axs[i,2].set_ylim(-.12,0.02)

fig.suptitle('Topoisomerases-Mediated DNA Relaxation and Supercoiling', fontsize=title_size)

# Collect
# exp_substrates.append(supercoiled_DNA)
# exp_products.append(relaxed_DNA)
# exp_superhelicals.append(sigma)


plt.savefig(file_out + '.png')
plt.savefig(file_out + '.pdf')
plt.show()
