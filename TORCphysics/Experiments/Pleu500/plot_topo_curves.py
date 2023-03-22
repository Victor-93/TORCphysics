import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from TORCphysics import EnvironmentFactory
from TORCphysics import effect_model as em

width = 6
height = 4

x = np.arange(-.1,.1,.001)
dt = 1
sites_filename = 'sites.csv'
environment_filename = 'environment.csv'
my_environment = EnvironmentFactory(environment_filename, [])
topo = my_environment.environment_list[0]
gyra = my_environment.environment_list[1]

topocurve = em.topo1_continuum(x, topo.concentration, topo.k_cat, dt)
gyracurve = em.gyrase_continuum(x, gyra.concentration, gyra.k_cat, dt)


fig, ax = plt.subplots(1, figsize=(width, height), tight_layout=True)

ax.plot(x,topocurve+gyracurve, 'black',  label='diff', lw=2)
ax.plot(x,topocurve, 'orange', label='topo I', lw=2)
ax.plot(x,gyracurve, 'purple', label='gyra', lw=2)

ax.legend(loc='best')
ax.grid(True)
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'$\sigma$ removed per timestep')
ax.set_title('Topoisomerase activity')
ax.set_xlim(x.min(),x.max())

plt.savefig("topo_curves.pdf")
plt.savefig("topo_curves.png")
