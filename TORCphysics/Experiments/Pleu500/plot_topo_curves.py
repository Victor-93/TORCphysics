import numpy as np
import matplotlib.pyplot as plt
from TORCphysics import EnvironmentFactory
from TORCphysics import effect_model as em
from TORCphysics.src import analysis as an

width = 6
height = 4

x = np.arange(-.1, .1, .001)
dt = 1
sites_filename = 'sites.csv'
environment_filename = 'environment.csv'
my_environment = EnvironmentFactory(environment_filename, [])
topo = my_environment.environment_list[0]
gyra = my_environment.environment_list[1]

topocurve = em.topo1_continuum(x, topo.concentration, topo.k_cat, dt)
gyracurve = em.gyrase_continuum(x, gyra.concentration, gyra.k_cat, dt)

fig, axs = plt.subplots(2, figsize=(width, 2 * height), tight_layout=True)

# Continuum
ax = axs[0]
ax.plot(x, topocurve + gyracurve, 'black', label='diff', lw=2)
ax.plot(x, topocurve, 'orange', label='topo I', lw=2)
ax.plot(x, gyracurve, 'purple', label='gyra', lw=2)
ax.legend(loc='best')
ax.grid(True)
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'$\sigma$ removed per timestep')
ax.set_title('Topoisomerase activity')
ax.set_xlim(x.min(), x.max())

# Stochastic
ax = axs[1]
topocurve, x = an.topoisomerase_activity_curves_stochastic(topo, sigma_min=-.2, sigma_max=.2)
gyracurve, x = an.topoisomerase_activity_curves_stochastic(gyra, sigma_min=-.2, sigma_max=.2)
# ax.plot(x, topocurve + gyracurve, 'black', label='diff', lw=2)
ax.plot(x, topocurve, 'orange', label='topo I', lw=2)
ax.plot(x, gyracurve, 'purple', label='gyra', lw=2)
ax.legend(loc='best')
ax.grid(True)
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'Rate')
ax.set_title('Topoisomerase activity stochastic')
# ax.set_xlim(x.min(), x.max())

plt.savefig("topo_curves.pdf")
plt.savefig("topo_curves.png")
