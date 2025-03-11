import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#sys_name = 'Ecoli'
#sys_name = 'Sal_noRNAPs'
#sys_name = 'no_topos'
sys_name='minimum_Sal'
#inputfile = 'torc_min_minimum_Ecoli_sites_df.csv'
#inputfile = 'torc_min_minimum_Sal_sites_df.csv'
#inputfile = 'torc_min_minimum_'+sys_name+'_sites_df.csv'
#enzymefile = 'torc_min_minimum_'+sys_name+'_enzymes_df.csv'
inputfile = 'torc_min_'+sys_name+'_sites_df.csv'
enzymefile = 'torc_min_'+sys_name+'_enzymes_df.csv'

df = pd.read_csv(inputfile)
edf = pd.read_csv(enzymefile)

tnames = ['PleuWT', 'tetA']#, 'lac1']
colors = ['yellow', 'purple', 'green']

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 20
height = 7*2
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16
# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we do the process
fig, axs = plt.subplots(2, figsize=(width, 1*height), tight_layout=True, sharex=True)

ax =axs[0]
ax.set_title(sys_name)

# Let's plot localsuperhelical
for i, name in enumerate(tnames):
    mask = df['name'] == name
    sigma = df[mask]['superhelical'].to_numpy()
    ax.plot(sigma, lw=1.5, color=colors[i], label=name)

# Let's plot binding/unbinding
for i, name in enumerate(tnames):
    mask = df['name'] == name
    binding_mask = df[mask]['binding'] == 1
    unbinding_mask = df[mask]['unbinding'] == 1
    if name == 'tetA':
        tdf = df[mask][binding_mask]
        ntdf = df[mask][unbinding_mask]
        a=2
    binding = df[mask]['binding'].to_numpy() - .5
    ax.plot(binding, lw=5, color=colors[i],alpha=0.3)


# Let's try to figure out elongation
for i, name in enumerate(tnames):
    mask = df['name'] == name
    time = df[mask]['time'].to_numpy()
#    binding = df[mask]['binding'].to_numpy()
#    enzymes = df[mask]['#enzymes'].to_numpy()
#    elongation = np.abs(binding-enzymes)-.5
    #axs.plot(time, elongation, lw=5, color=colors[i],alpha=0.3)
#    axs.fill_between(time,np.ones_like(elongation)*(-.5), elongation, lw=1, color=colors[i],alpha=0.3)
circuit_mask = df['type'] == 'circuit'

# Let's try to figure out elongation but from the enzymes file
for i, name in enumerate(tnames):
    new_df = df[circuit_mask]

    new_df['elongation'] = 0

    mask = edf['site'] == name
    objects = edf[mask]
    mask_elongation = objects['name'] == 'RNAP_Elongation'
    objects = objects[mask_elongation]
    objects.drop_duplicates(subset=['time'],inplace=True)

    elongation_frames = objects['time'].unique()  # Assuming 'time' represents frame values
    for f in elongation_frames:
        p = int(f)
        new_df['elongation'].iloc[p] = 1
    elongation = new_df['elongation'].to_numpy() - .5
#    new_df.loc[new_df['frame'].isin(elongation_frames), 'elongation'] = 1
    ax.fill_between(time,np.ones_like(elongation)*(-.5), elongation, lw=1, color=colors[i],alpha=0.4)
    axs[1].fill_between(time,np.ones_like(elongation)*(-.5), elongation*50, lw=1, color=colors[i],alpha=.4)

# Let's put some other with very shaded areas where the RNAPs bind but there's no elongation
for i, name in enumerate(tnames):
    new_df = df[circuit_mask]
    new_df['elongation'] = 0

    mask = edf['site'] == name
    objects = edf[mask]
    mask_elongation = objects['name'] == 'RNAP'
    objects = objects[mask_elongation]
    objects.drop_duplicates(subset=['time'],inplace=True)

    elongation_frames = objects['time'].unique()  # Assuming 'time' represents frame values
    for f in elongation_frames:
        p = int(f)
        new_df['elongation'].iloc[p] = 1
    elongation = new_df['elongation'].to_numpy() - .5
#    new_df.loc[new_df['frame'].isin(elongation_frames), 'elongation'] = 1
    ax.fill_between(time,np.ones_like(elongation)*(-.5), elongation, lw=1, color=colors[i],alpha=0.2)


# Let's plot the global
mask = df['type'] == 'circuit'
sigma = df[mask]['superhelical'].to_numpy()
ax.plot(sigma, lw=1.75,color='black', label='global')

ax.set_ylim([-0.12,0.0])
ax.set_xlabel('Time')
ax.set_ylabel('Local superhelical')
ax.legend(loc='best')
ax.grid(True)

#plt.savefig('supercoiling.png')
#plt.savefig(out_file+'png')

# Number of topos
# ----------------------------------------------------------------------------------------------------------------------
ax =axs[1]

topo_names = ['DNA_topoI', 'DNA_gyrase']
topo_colors = ['red', 'cyan']
for i, name in enumerate(topo_names):
    mask = df['name'] == name
    time = df[mask]['time'].to_numpy()
    enzymes = df[mask]['#enzymes'].to_numpy()
    ax.plot(time, enzymes, lw=1, color=topo_colors[i], label=name)

ax.set_xlabel('Time')
ax.set_ylabel('Number of Enzymes')
ax.legend(loc='best')
ax.grid(True)
ax.set_ylim([0,15])


plt.show()




