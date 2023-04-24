import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
import pandas as pd
import params
import statistical_model as sm
import sys


width = 6
height = 4

globalsigma = np.loadtxt('../../global_sigma.txt')
t = globalsigma[:,1]
s = globalsigma[:,0]

fig, ax = plt.subplots(1, figsize=(width, height), tight_layout=True)

ax.plot(t,s, 'red')

#ax.legend(loc='best')
ax.grid(True)
ax.set_ylabel(r'$\sigma$')
ax.set_xlabel(r'time (seconds)')

plt.savefig("sigma_curve.pdf")
