import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
from TORCphysics import utils
from TORCphysics import binding_model as bm

# Description
#-----------------------------------------------------------------------------------------------------------------------
# Plots detailed parameters for a particular case. It includes susceptibility, prod rate, number of enzymes and
# local/global superhelical density

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
#promoter_cases = ['weak', 'medium', 'strong']
promoter_cases = ['medium']
dt=1.0#0.5

# Selec this for V2
path = 'genearch_V2/'
model_code='sus_genearch_V2-01-'
figtitle = model_code + promoter_cases[0]
outfile = model_code+'trials'

# Select this one for V1
#path = 'genearch_V1/'
#model_code='sus_genearch_V1-01-'
#figtitle = model_code + promoter_cases[0] + ' gamma0.1'
#outfile = model_code+'trials'

# Select this one for V0
#path = 'genearch_V0/'
#model_code='sus_genearch_V0-01-'
#figtitle = model_code + promoter_cases[0] + ' gamma0.05'
#outfile = model_code+'trials'

susceptibility_files = []
trials_file = []
for pcase in promoter_cases:
    susceptibility_files.append('../junier_data/'+ pcase + '.csv')

    # This is for V2
    trials_file.append(path+model_code+pcase+'_dt'+str(dt)+'-trials.pkl')

    # This for V1 with gamma
    # trials_file.append(path+model_code+pcase+'_gamma0.1_dt'+str(dt)+'-trials.pkl')

    # This for V0 with gamma
    #trials_file.append(path+model_code+pcase+'_gamma0.05_dt'+str(dt)+'-trials.pkl')

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 7
height = 3.75
lw = 4
font_size = 12
xlabel_size = 14
title_size = 16

# line styles
model_ls = '-o'
exp_ls = '--o'
titles = ['Weak Promoter', 'Medium Promoter', 'Strong Promoter']

colors = ['green', 'blue', 'red']

# Processing functions
#-----------------------------------------------------------------------------------------------------------------------
# Load
#-----------------------------------------------------------------------------------------------------------------------
pickle_data = []
for pickle_file in trials_file:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_data.append(data)

# Plot
#-----------------------------------------------------------------------------------------------------------------------
#sigma = np.arange(-.2, .0, 0.001)
# Let's plot as we do the process
fig, axs = plt.subplots(3,2, figsize=(2*width, 3*height), tight_layout=True, sharex=True)

fig.suptitle(figtitle, fontsize=title_size)

outside_label = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)']

color_on = 'black'
color_off = 'red'

color_enzyme = {'RNAP': 'black', 'topoI':'red', 'gyrase':'blue'}
nenzyme_ylim = [0,6]

n = len(trials_file)
for i in range(n):

    # Load data
    # ----------------------------------------------------------------------------------
    data_list = pickle_data[i].results

    # Filter data_list to include only entries with 'status' == 'ok'
    data_list = [item for item in data_list if item['status'] == 'ok']

    # Find the dictionary with the smallest 'loss'
    data = min(data_list, key=lambda x: x['loss'])
    x = data['distance']

    # Susceptibility
    # ----------------------------------------------------------------------------------
    ax = axs[0,0]
    pr = data['prod_rate'][:,0]
    prs = data['prod_rate'][:,1]/np.sqrt(60) # Error in the mean
    y = pr/pr[4] #susceptibility
    ys = np.zeros_like(y)
    for j in range(len(y)):
        print(j)
        ys[j] = y[j] * np.sqrt( np.square(prs[j]/pr[j]) + np.square(prs[4]/pr[4]) )
#    y = data['susceptibility'][:,0]
#    ys = data['susceptibility'][:,1]

    ax.set_title('Susceptibility ', fontsize=title_size)

    ax.plot(x, y, model_ls, lw=lw, color='black')
    ax.fill_between(x, y-ys, y+ys, lw=lw, color='black', alpha=0.2)

    ax.set_ylabel(r'Susceptiblity', fontsize=font_size)
    ax.grid(True)
    ax.set_ylim([.7, 1.6])


    # Production rate
    # ----------------------------------------------------------------------------------
    ax = axs[0,1]
    y = data['prod_rate'][:,0]
    ys = data['prod_rate'][:,1]

    ax.set_title('Production Rate', fontsize=title_size) #color=colors[i],

    ax.plot(x, y, model_ls, lw=lw, color='black')
    ax.fill_between(x, y-ys, y+ys, lw=lw, color='black', alpha=0.2)

    ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    ax.set_ylabel(r'Production rate $(s^{-1})$', fontsize=font_size)
    ax.grid(True)
    ax.set_xlim([80, 6000])
    ax.set_xscale('log')

    # Number of bound Topoisomerase I
    # ----------------------------------------------------------------------------------
    ax = axs[1,0]
    ax.set_title('Number of Bound Topoisomerases', fontsize=title_size) #color=colors[i],
    colors = ['red', 'blue']
    for j, name in enumerate(['topoI', 'gyrase']):
        y = data['nenzymes'][name][:,0]
        ys = data['nenzymes'][name][:,1]

        ax.errorbar(x,y,yerr=ys, fmt='-o', color=colors[j], label=name)

    ax.set_ylabel(r'Bound Molecules', fontsize=font_size)
    ax.grid(True)
    #ax.set_ylim(nenzyme_ylim)

    # Number of bound RNAPs
    # ----------------------------------------------------------------------------------
    ax = axs[1,1]
    name = 'RNAP'
    y = data['nenzymes'][name][:,0]
    ys = data['nenzymes'][name][:,1]


    ax.set_title('Number of Bound RNAPs', fontsize=title_size) #color=colors[i],

    ax.errorbar(x, y, yerr=ys, fmt='-o', color='black', label=name)

    ax.set_ylabel(r'Bound Molecules', fontsize=font_size)
    ax.grid(True)

    # Local superhelical
    # ----------------------------------------------------------------------------------
    ax = axs[2,0]
    y = data['local_superhelical'][:,0]
    ys = data['local_superhelical'][:,1]

    ax.set_title('Superhelicity at the Promoter', fontsize=title_size) #color=colors[i],

    ax.errorbar(x,y,yerr=ys, fmt='-o', color='black')

    ax.grid(True)
    ax.set_ylabel('Superhelical Density', fontsize=font_size)
    ax.set_ylim([-.1,-0.01])

    # Global superhelical
    # ----------------------------------------------------------------------------------
    ax = axs[2,1]
    y = data['global_superhelical'][:,0]
    ys = data['global_superhelical'][:,1]

    ax.set_title('Global Superhelicity', fontsize=title_size) #color=colors[i],

    ax.errorbar(x,y,yerr=ys, fmt='-o', color='black')

    ax.set_ylabel('Superhelical Density', fontsize=font_size)
    ax.set_ylim([-.1,.05])

    ax.set_xlabel('Upstream Distance (bp)', fontsize=font_size)
    ax.grid(True)
    ax.set_xlim([80, 6000])
    ax.set_xscale('log')


ax_list = [axs[0,0], axs[0,1], axs[1,0], axs[1,1], axs[2,0], axs[2,1]]
for i, ax in enumerate(ax_list):
    ax.legend(loc='upper left')
#    ax.legend(loc='best')

    # Add label outside the plot
    ax.text(-0.1, 1., outside_label[i], transform=ax.transAxes,
            fontsize=font_size*1.5, fontweight='bold', va='center', ha='center')

#plt.savefig(outfile+'.pdf')
plt.show()
