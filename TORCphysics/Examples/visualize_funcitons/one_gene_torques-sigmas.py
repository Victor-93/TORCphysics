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
fig, axs = plt.subplots(4, figsize=(width, 4 * height), tight_layout=True)

# First for sigma
x0 = 300
lg = 900
xd = 320
sigma0 = -0.046
v0 = 30.0
gamma = .1

total_time = lg/v0
#total_time =
time = np.arange(0,total_time,0.1)
sigmaup = sigma_upstream_v0(sigma0, x0, v0, gamma, time)
sigmadown=sigma_downstream_v0(sigma0,xd,lg,v0,gamma,time)

# sigma
ax = axs[0]
ax.plot(time, sigmaup, lw=lw, color=up_c,label='upstream')
ax.plot(time, sigmadown, lw=lw, color=down_c, label='downstream')
ax.grid(True)
ax.legend(loc='best')
ax.set_xlabel('Time (s)')
ax.set_ylabel(r'$\sigma$')

# Let's plot the torque
# Let's just define a site_list to use for various tests
site_gene1 = Site(site_type='gene', name='test_gene1', start=100, end=500, k_on=3.00)
site_gene2 = Site(site_type='gene', name='test_gene2', start=600, end=800, k_on=3.00)
site_gene3 = Site(site_type='gene', name='test_gene3', start=1200, end=1000, k_on=3.00)
site_tf = Site(site_type='TF', name='test_TF', start=1200, end=1000, k_on=3.00)
site_list1 = [site_gene1, site_gene2, site_gene3, site_tf]

# And list of enzymes (RNAPs without effect model)
RNAP = Enzyme(e_type='RNAP', name='test1', site=site_list1[0], size=30, effective_size=15, position=30,
              twist=0.0, superhelical=0.0, effect_model_name='RNAPStall')

ax = axs[1]
torque_up = np.zeros_like(time)
torque_down = np.zeros_like(time)
for i,s in enumerate(sigmaup):
    torque_up[i] = ef.Marko_torque(s)
for i,s in enumerate(sigmadown):
    torque_down[i] = ef.Marko_torque(s)
torque_T = torque_down - torque_up
ax.plot(time, torque_up, lw=lw, color=up_c,label='upstream')
ax.plot(time, torque_down, lw=lw, color=down_c, label='downstream')
ax.plot(time, torque_T, lw=lw, color='black', label='total')
ax.grid(True)
ax.legend(loc='best')
ax.set_xlabel('Time (s)')
ax.set_ylabel(r'Torque (pN)')

# Let's do a heatmap
# ---------------------------------------
ax = axs[2]

# Torque
sigma_up = np.arange(-.15,.151,0.001)
sigma_down = np.arange(-.15,.151,0.001)

torque = np.zeros( (len(sigma_up), len(sigma_down)))
for i, s_up  in enumerate(sigma_up):
    for j, s_down in enumerate(sigma_down):
        t_up = ef.Marko_torque(s_up)
        t_down = ef.Marko_torque(s_down)
        total_torque = t_down - t_up
        #if abs(total_torque) >= params.stall_torque:
        if total_torque >= params.stall_torque:
                torque[i,j] = 0
        else:
            torque[i,j] = 1
        #torque[i,j]= total_torque

print('sigma_s', params.sigma_s, ' sigma_p', params.sigma_p)

# Display the heatmap
cax = ax.imshow(torque, aspect='auto', cmap='Greens_r', origin='lower', extent=[sigma_down.min(), sigma_down.max(), sigma_up.min(), sigma_up.max()])

# Add a colorbar
cbar = fig.colorbar(cax)
cbar.set_label('0-stall, 1-moves')

# Set axis labels
ax.set_xlabel('sigma_down')
ax.set_ylabel('sigma_up')

# Plot the evolution of our simulated sigmaup and down
ax.plot(sigmadown,sigmaup)

# Let's draw the sigma_p (the regime where plectonemes dominate and RNAP stalls)
print('sigma_s', params.sigma_s, ' sigma_p', params.sigma_p)
ax.plot([-params.sigma_p, params.sigma_p],[params.sigma_p, params.sigma_p], 'red')
ax.plot([-params.sigma_p, params.sigma_p],[-params.sigma_p, -params.sigma_p], 'red')
ax.plot([-params.sigma_p, -params.sigma_p],[-params.sigma_p, params.sigma_p], 'red')
ax.plot([params.sigma_p, params.sigma_p],[-params.sigma_p, params.sigma_p], 'red')

# Set x and y ticks to show actual values
ax.set_xticks(np.linspace(sigma_down.min(), sigma_down.max(), num=5))
ax.set_yticks(np.linspace(sigma_up.min(), sigma_up.max(), num=5))

ax.grid(True)

# Let's try to do an experiment
# ---------------------------------------
ax = axs[3]
dx = 100
x = 300 # initial x
iterations = 50
lg = 900 #let's increase the gene length...
total_time = lg/v0
time = np.arange(0,total_time,0.1)

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

if len(t_list) == 0:
    t_list = np.zeros_like(x_list)

ax.plot(x_list, t_list)

ax.set_xscale('log')
ax.grid(True)
ax.set_xlabel('upstream distance (bp)')
ax.set_ylabel('stall time')

plt.show()
