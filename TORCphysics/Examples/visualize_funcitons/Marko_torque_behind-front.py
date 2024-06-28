from TORCphysics import effect_model as ef
from TORCphysics import Enzyme, Site, params
import matplotlib.pyplot as plt
import numpy as np

# What I want to do here, is to visualize Marko's elasticity model, but, this time, we assue two torques that act
# on an enzyme (RNAP in this case). We have a torque behind (upstream) and torque ahead (downstream).

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

# Plot
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, figsize=(width, 1*height), tight_layout=True)

# Torque
sigma_up = np.arange(-.2,.2,0.001)
sigma_down = np.arange(-.2,.2,0.001)

torque = np.zeros( (len(sigma_up), len(sigma_down)))
for i, s_up  in enumerate(sigma_up):
    for j, s_down in enumerate(sigma_down):
        t_up = ef.Marko_torque(s_up)
        t_down = ef.Marko_torque(s_down)
        total_torque = t_down - t_up
        if total_torque >= params.stall_torque:
            torque[i,j] = 0
        else:
            torque[i,j] = 1

print('sigma_s', params.sigma_s, ' sigma_p', params.sigma_p)

ax = axs
# Display the heatmap
cax = ax.imshow(torque, aspect='auto', cmap='Greens_r', origin='lower', extent=[sigma_down.min(), sigma_down.max(), sigma_up.min(), sigma_up.max()])

# Add a colorbar
cbar = fig.colorbar(cax)
cbar.set_label('0-stall, 1-moves')

# Set axis labels
ax.set_xlabel('sigma_down')
ax.set_ylabel('sigma_up')

# Set x and y ticks to show actual values
ax.set_xticks(np.linspace(sigma_down.min(), sigma_down.max(), num=5))
ax.set_yticks(np.linspace(sigma_up.min(), sigma_up.max(), num=5))

# Show the plot
plt.show()