from TORCphysics import effect_model as ef
from TORCphysics import Enzyme, Site, params, utils
import matplotlib.pyplot as plt
import numpy as np

# What I want to do here, is to visualize Marko's elasticity model, and plot it's torque

# Plotting params
# ----------------------------------------------------------------------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16
sigma_s_c = 'red'
sigma_p_c = 'blue'
torque_color = 'black'

C_length = params.C_length
A_length = params.A_length
f_stretching = params.f_stretching
kBT_pN_nm = params.kBT_pN_nm

# Print values of interes
print('p',params.p_stiffness)
print('c', params.c_stiffness)
print('cs',params.cs_energy)
print(1.0 - ((C_length / (4.0 * A_length)) *
             np.sqrt(kBT_pN_nm / (A_length * f_stretching))))
print(params.g_energy)
print(params.P_length)

# Plot
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(2, figsize=(width, 2*height), tight_layout=True)

# Torque
sigma = np.arange(-.2,.2,0.0001)
torque = np.zeros(len(sigma))
print('sigma_s', params.sigma_s, ' sigma_p', params.sigma_p)

for i,s in enumerate(sigma):
    torque[i] = utils.Marko_torque(s)

ax = axs[0]
ax.plot(sigma, torque, lw=lw, color=torque_color)
ax.plot([params.sigma_s,params.sigma_s], [-100,100], lw=lw, color=sigma_s_c)
ax.plot([-params.sigma_s,-params.sigma_s], [-100,100], lw=lw, color=sigma_s_c)
ax.plot([params.sigma_p,params.sigma_p], [-200,200], lw=lw, color=sigma_p_c)
ax.plot([-params.sigma_p,-params.sigma_p], [-200,200], lw=lw, color=sigma_p_c)
#ax.set_ylim(-40,40)
ax.set_xlim(-.1,.1)
ax.grid(True)
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'torque (pN)')

# Let's plot the velocity
# Let's just define a site_list to use for various tests
site_gene1 = Site(site_type='gene', name='test_gene1', start=100, end=500, k_on=3.00)
site_gene2 = Site(site_type='gene', name='test_gene2', start=600, end=800, k_on=3.00)
site_gene3 = Site(site_type='gene', name='test_gene3', start=1200, end=1000, k_on=3.00)
site_tf = Site(site_type='TF', name='test_TF', start=1200, end=1000, k_on=3.00)
site_list1 = [site_gene1, site_gene2, site_gene3, site_tf]

# And list of enzymes (RNAPs without effect model)
RNAP = Enzyme(e_type='RNAP', name='test1', site=site_list1[0], size=30, effective_size=15, position=30,
                 twist=0.0, superhelical=0.0, effect_model_name='RNAPStall')

torque2= np.arange(-50,50,.1)
v = ef.velocity_2022SevierBioJ(RNAP, torque)
v2 = ef.velocity_2022SevierBioJ(RNAP, torque2)
ax = axs[1]
ax.plot(torque, v, '-o', color='blue')
ax.plot(torque2, v2, '-o', color=torque_color)
ax.grid(True)
#ax.set_ylim(0,80)
#ax.set_xlim(-100,100)

plt.show()