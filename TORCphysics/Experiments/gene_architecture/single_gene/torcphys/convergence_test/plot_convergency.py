import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from TORCphysics import analysis as ann

#TIME TO PLOT

# **********************************************************************************************************************
# Inputs/Initial conditions
# **********************************************************************************************************************
n_simulations = 4 # 16

idist = 0#11 #Which distance we want
promoter_case = 'weak'

# Junier experimental data - we only need distances for now.
# --------------------------------------------------------------
experimental_file = '../../junier_data/' + promoter_case + '.csv'

# Promoter responses
# --------------------------------------------------------------
promoter_response = '../../promoter_responses/'+promoter_case+'.csv'

# Simulation conditions
# --------------------------------------------------------------

# The system is like:  UP______GENE____DOWN, where UP and DOWN are upstream and downstream barriers.
# The upstream barrier will always be located at 0 (linear system), the gene start site is located at
# x=upstream_distance and termination at upstream_distance+gene_length, and the total length of the region is
# upstream_distance + gene_length + downstream_distance
gene_length = 900
downstream_distance = 320

# Initial superhelical density
sigma0 = -0.046

dt = 0.25
initial_time = 0
final_time = 20000#36000#3600 #9000 ~2.5hrs
time = np.arange(initial_time, final_time + 2*dt, dt)
frames = len(time)


# Circuit conditions
# --------------------------------------------------------------
output_prefix = 'single_gene'  #Although we won't use it
series = True
continuation = False
enzymes_filename = None
environment_filename = 'environment.csv'
# We will define sites and circuit csv as we advance

# These two filenames will be created and updated for every single cases
circuit_filename = promoter_case+'-circuit.csv'
sites_filename = promoter_case+'-sites.csv'

# Plotting params
# ----------------------------------------------------------------------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16
exp_color = 'green'
model_color = 'blue'
torc_color = 'red'
weak_color = 'green'
medium_color = 'blue'
strong_color = 'red'
ms = 10

# **********************************************************************************************************************
# Load data and prepare systems
# **********************************************************************************************************************
# Load exp
experimental_curve = pd.read_csv(experimental_file)

# Extract distances
distances = list(experimental_curve['distance'])

# **********************************************************************************************************************
# Functions
# **********************************************************************************************************************
#def cal_steady(y):
#
#    frames = len(y)
##    accumulation = np.zeros(frames)  # this is the sum of the events
 #   log_accumulation = np.zeros(frames)  # And the logarithm of these sum of events divided by the time.
 #   for k in range(1, frames):
 #       accumulation[k] = np.sum(y[0:k])
##
 #   if np.sum(accumulation) <= 0.0:
 #       rate = 0.0000000001
 #   else:
  ##      log_accumulation = np.log(accumulation / time)
  #      plateau = np.mean(log_sum_ini[ta:tb])
  #      rate = np.exp(plateau)
  ##  curve = log_sum_ini
ylabels = ['binding', 'unbinding', r'global $\sigma$', r'local $\sigma$', 'kde']
enzyme_color = {'RNAP':'black', 'topoI':'red', 'gyrase':'cyan'}
# Plot
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(5, figsize=(width, 5*height), tight_layout=True)

for i, upstream_distance in enumerate(distances):
    if i != idist:
        continue

    fig.suptitle('distance '+ str(upstream_distance))
    gbinding = np.zeros_like(time)
    gunbinding = np.zeros_like(time)
    gsigma = np.zeros_like(time)
    superhelical = np.zeros_like(time)

    binding_df = pd.read_csv(promoter_case+'_binding_distance-' + str(int(upstream_distance))+'.csv')
    unbinding_df = pd.read_csv(promoter_case+'_unbinding_distance-' + str(int(upstream_distance))+'.csv')
    superhelical_df = pd.read_csv(promoter_case+'_superhelical_distance-' + str(int(upstream_distance))+'.csv')
    gsigma_df = pd.read_csv(promoter_case+'_gsigma_distance-' + str(int(upstream_distance))+'.csv')

    # Binding
    ax = axs[0]
    for column_name, column_data in binding_df.items():
        y = column_data.to_numpy()  # Convert column data to a NumPy array
        gbinding = gbinding + y
        sum_time, log_sum, rate = ann.calculate_steady_state_initiation_curve(y, time, ta=int(frames/2), tb=-1)
        ax.plot(time,sum_time)
        print('bind',rate)
    y = gbinding/n_simulations
    sum_time, log_sum, rate_bind = ann.calculate_steady_state_initiation_curve(y, time, ta=int(frames / 2), tb=-1)
    ax.plot(time, sum_time, '-k')

    # Unbinding
    ax = axs[1]
    for column_name, column_data in unbinding_df.items():
        y = column_data.to_numpy()  # Convert column data to a NumPy array
        gunbinding = gunbinding + y
        sum_time, log_sum, rate = ann.calculate_steady_state_initiation_curve(y, time, ta=int(frames/2), tb=-1)
        ax.plot(time,sum_time)
        print('unbind',rate)
    y = gunbinding/n_simulations
    sum_time, log_sum, rate_prod = ann.calculate_steady_state_initiation_curve(y, time, ta=int(frames / 2), tb=-1)
    ax.plot(time, sum_time, '-k')

    # global supercoilin
    ax = axs[2]
    for column_name, column_data in gsigma_df.items():
        y = column_data.to_numpy()  # Convert column data to a NumPy array
        gsigma = gsigma + y
        #ax.plot(time,y)
    gsigma = gsigma/n_simulations
    mean = np.mean(gsigma)
    ax.plot(time, gsigma, '-r')
    ax.plot([time[0], time[-1]], [mean,mean], '-k')

    # supercoilint at site
    ax = axs[3]
    for column_name, column_data in superhelical_df.items():
        y = column_data.to_numpy()  # Convert column data to a NumPy array
        superhelical = superhelical + y
        #ax.plot(time,y)
    superhelical = superhelical/n_simulations
    mean = np.mean(superhelical)
    ax.plot(time, superhelical, '-r')
    ax.plot([time[0], time[-1]], [mean,mean], '-k')


    # KDEs
    ax = axs[4]
    for name in ['RNAP', 'topoI', 'gyrase']:
        kde = np.loadtxt(promoter_case+'_'+name+'_KDE_distance-'+str(int(upstream_distance))+'.csv')
        ax.plot(kde[:,0], kde[:,1]/np.max(kde[:,1]), '-', color=enzyme_color[name])


    print('rate_binding', rate_bind, 'rate_prod', rate_prod)
    break

#y = gbinding
#sum_time, log_sum, rate = ann.calculate_steady_state_initiation_curve(y, time, ta=int(frames / 2), tb=-1)
#axs[0].plot(time,sum_time, '-k')
#y = gunbinding
#sum_time, log_sum, rate = ann.calculate_steady_state_initiation_curve(y, time, ta=int(frames / 2), tb=-1)
#axs[1].plot(time,sum_time, '-k')

for n in range(5):
    axs[n].grid(True)
    axs[n].set_xlabel('Time (s)')
    axs[n].set_ylabel(ylabels[n], fontsize=xlabel_size)
#axs[0].set_ylim(0,0.01)
axs[3].set_ylim(-0.2,0.0)
axs[4].set_xlabel('region (bp)')

# Define the text to display
trate_bind = f'kbind={rate_bind:.4f}'
trate_prod = f'kprod={rate_prod:.4f}'

# Add the text box to the plot
props = dict(boxstyle='round', facecolor='silver', alpha=0.2)
axs[0].text(0.5, 0.8, trate_bind, transform=axs[0].transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
axs[1].text(0.5, 0.1, trate_prod, transform=axs[1].transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

#axs[1].set_ylabel('Medium promoter', fontsize=xlabel_size)
#a#xs[2].set_ylabel('Strong promoter', fontsize=xlabel_size)
#axs[2].set_xlabel('Distance (bp)', fontsize=xlabel_size)
#a#xs[2].legend(loc='best', fontsize=font_size)

plt.savefig(promoter_case+'-dist-'+str(idist)+'.png')
#plt.savefig('susceptibility.png')
#plt.savefig('junier.pdf')
plt.show()

