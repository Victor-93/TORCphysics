from TORCphysics import effect_model as ef
from TORCphysics import Enzyme, Site, params
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm


# So, we have a system with a single gene and a RNAP transcribing.
# Here, we want to see how the torques, velocities and sigmas will evolve as the RNAP advances.

#  UP-------->>>>>>>>>>>>>>>---------DOWN
# We have a UPstream barrier, a gene, and then a downstream barrier.


# sigma0 = initial superhelical density
# x0 = distance from transcription start site and upstream barrier
# v0 = RNAP velocity
# gamma = supercoiling injection
# xd = distance between termination site and downstream barrier
# lg = gene length
def sigma_upstream_v0(sigma0, x0, v0, gamma, time):
    sigma = np.zeros_like(time)
    for i, t in enumerate(time):
        sigma[i] = (sigma0 * x0 - gamma * v0 * t) / (x0 + v0 * t)
    return sigma


def sigma_downstream_v0(sigma0, xd, lg, v0, gamma, time):
    sigma = np.zeros_like(time)
    for i, t in enumerate(time):
        sigma[i] = (sigma0 * (lg + xd) + gamma * v0 * t) / (lg + xd - v0 * t)
    return sigma


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
up_c = 'red'
down_c = 'blue'

# Plot
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, figsize=(width, 1 * height), tight_layout=True)

# First for sigma
x0 = 300
lg = 900
xd = 320
sigma0 = -0.046
v0 = 30.0

# Let's try to do an experiment
# ---------------------------------------
ax = axs
dx = 100
x0 = 100
iterations = 50
lg = 900 #let's increase the gene length...
total_time = lg/v0
time = np.arange(0,total_time,0.1)

gammas = np.arange(0.1, 1.001, 0.001)

# Choose a colormap
colormap = cm.inferno

# Normalize the gamma values to the range [0, 1] for the colormap
norm = plt.Normalize(gammas.min(), gammas.max())

# Plot each line with a color from the colormap
for gamma in gammas:

    x = x0  # initial x
    torque_up = np.zeros_like(time)
    torque_down = np.zeros_like(time)
    x_list = []
    t_list = []

    for j in range(iterations):
        x_list.append(x)

        sigmaup = sigma_upstream_v0(sigma0, x, v0, gamma, time)
        sigmadown = sigma_downstream_v0(sigma0, xd, lg, v0, gamma, time)

        for i, s in enumerate(sigmaup):
            torque_up[i] = ef.Marko_torque(s)
        for i, s in enumerate(sigmadown):
            torque_down[i] = ef.Marko_torque(s)

        total_torque = torque_down - torque_up

        for k, t in enumerate(total_torque):
            if t >= params.stall_torque:
                t_list.append(time[k])
                break
        x = x + dx

    if len(t_list) != len(x_list):
        t_list = np.zeros_like(x_list)
        print(gamma)

    color = colormap(norm(gamma))
    ax.plot(x_list, t_list, color=color)


# Create a scalar mappable for the colorbar
sm = cm.ScalarMappable(cmap=colormap, norm=norm)
sm.set_array([])  # Only needed for matplotlib < 3.1

# Add the colorbar to the plot
cbar = fig.colorbar(sm, ax=ax)
cbar.set_label('Gamma')

ax.set_xscale('log')
ax.grid(True)
ax.set_xlabel('upstream distance (bp)')
ax.set_ylabel('stall time')

plt.show()
