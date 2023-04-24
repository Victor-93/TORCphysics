import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
import pandas as pd
import params
import statistical_model as sm
import sys

#topo_w = params.topo_w
#topo_t = params.topo_t
#topo_c = params.topo_c
#topo_k = params.topo_k

#gyra_w = params.gyra_w
#gyra_t = params.gyra_t
#gyra_c = params.gyra_c
#gyra_k = params.gyra_k
width = 6
height = 4

x = np.arange(-.15,.15,.001)

topo = sm.topoisomerase_activity(x)
gyra = sm.gyrase_activity(x)

#topo = np.abs(topo)
#topo = topo/np.max(topo)
#gyra = np.abs(gyra)
#gyra = gyra/np.max(gyra)

fig, ax = plt.subplots(1, figsize=(width, height), tight_layout=True)

#ax.plot(x,topo+gyra, 'black',  label='diff')
#ax.plot(x,topo, 'orange', label='topo I')
#ax.plot(x,gyra, 'purple', label='gyra')
#ax.plot(x,topo2+gyra2, 'black',  label='diff')
#ax.plot(x,topo2, 'orange', label='topo I')
#ax.plot(x,gyra2, 'purple', label='gyra')

ax.plot(x,topo+gyra, 'black',  label='diff', lw=2)
ax.plot(x,topo, 'orange', label='topo I', lw=2)
ax.plot(x,gyra, 'purple', label='gyra', lw=2)

ax.legend(loc='best')
ax.grid(True)
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'$\sigma$ removed per timestep')
ax.set_title('Topoisomerase activity')
ax.set_xlim(-.1,.1)

plt.savefig("topo_curves.pdf")
plt.savefig("topo_curves.png")
