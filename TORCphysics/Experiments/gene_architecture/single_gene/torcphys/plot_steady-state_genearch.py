import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
from TORCphysics import utils
from TORCphysics import binding_model as bm

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
dt=1.0#0.5

k_weak=0.003
# k_weak=0.0334

model_code='reproduce-sus_GB-Stages-avgx2-01-'
model_code='sus_GB-Stages-avgx2-02-'

distance_list = [0, 1,2,3,4,5,6,7,8]  #List of distances to plot
ndist=len(distance_list)

experimental_files = []
susceptibility_files = []
calibration_files = []
params_files = []
for pcase in promoter_cases:
    experimental_files.append('../junier_data/inferred-rate_kw'+str(k_weak)+'_' + pcase + '.csv')
    susceptibility_files.append('../junier_data/'+ pcase + '.csv')

    # Susceptibility
    calibration_files.append('susceptibility/'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    #calibration_files.append('susceptibility/reproduce-'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    #calibration_files.append('susceptibility/reproduce-'+model_code+pcase+'_dt'+str(dt)+'_RNAPTracking_off.pkl')
    #calibration_files.append('susceptibility/reproduce-gene-avg_'+model_code+pcase+'_dt'+str(dt)+'.pkl')
    params_files.append('susceptibility/'+model_code+pcase+'_dt'+str(dt)+'.csv')

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 7
height = 3.5
lw = 1
font_size = 12
xlabel_size = 14
title_size = 16

# line styles
titles = ['Weak Promoter', 'Medium Promoter', 'Strong Promoter']
colors = {'RNAP': 'black', 'topoI': 'red', 'gyrase': 'blue'}

# Processing functions
#-----------------------------------------------------------------------------------------------------------------------
# Calculate rates as a function of distance
def get_nenzymes(results_list, dist_index):
    x = results_list['distance'][dist_index]
    nRNAPs = results_list['nenzymes']['RNAP'][dist_index]
    ntopoIs = results_list['nenzymes']['topoI'][dist_index]
    ngyrases = results_list['nenzymes']['gyrase'][dist_index]
    return {'RNAP': nRNAPs, 'topoI': ntopoIs, 'gyrase': ngyrases, 'distance': x}

# Load
#-----------------------------------------------------------------------------------------------------------------------
pickle_data = []
for pickle_file in calibration_files:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_data.append(data)

# Plot
#-----------------------------------------------------------------------------------------------------------------------
sigma = np.arange(-.13, .0, 0.001)
#sigma = np.arange(-.2, .0, 0.001)
# Let's plot as we do the process
fig, axs = plt.subplots(ndist,3, figsize=(2*width, ndist*height), tight_layout=True, sharex=False)

outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

for i, dist in enumerate(distance_list):

    for j, data in enumerate(pickle_data):
        ax = axs[i, j]

        e_dict = get_nenzymes(data[0], dist)
        distance = e_dict['distance']
        frames = len(e_dict['RNAP'][0])
        time = np.arange(0, frames)

        ax.set_title(titles[j] + ' upstream: ' + str(distance), fontsize=title_size) #color=colors[i],

        if i == 0:
            for name in ['gyrase', 'topoI', 'RNAP']:
                y = e_dict[name][0]
                ys = e_dict[name][1]
                ax.plot(time, y, lw=lw, color=colors[name], label=name)
        else:
            for name in ['gyrase', 'topoI', 'RNAP']:
                y = e_dict[name][0]
                ys = e_dict[name][1]
                ax.plot(time, y, lw=lw, color=colors[name])

        for name in ['gyrase', 'topoI', 'RNAP']:
            y = e_dict[name][0]
            ys = e_dict[name][1]
            ax.fill_between(time, y-ys, y+ys, lw=lw, color=colors[name], alpha=0.2)

        ax.set_xlabel('Time (seconds)', fontsize=font_size)
        ax.set_ylabel(r'Number of enzymes', fontsize=font_size)
        ax.grid(True)
        #ax.set_xlim([80, 6000])
        #ax.set_xscale('log')

axs[0,0].legend(loc='best')


plt.show()



