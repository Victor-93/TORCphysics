import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Description
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot the losses, and calculate averages. Only for the Stages-stall

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
percentage_threshold = .1
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0 #0.25
file_out = 'loss-small-dist_'+str(dt)
#param_out = 'avg_' # prefix for outputting averaged parameters to be used as model parametrization
param_out = 'avgx2_' # prefix for outputting averaged parameters to be used as model parametrization

#path = 'track-StagesStall/'
#path = 'track-StagesStall/avg02/'
path = '../calibration/'
loss_file = '-calibration_avg-RNAPTracking_nsets_p2_small_dt'+str(dt)+'-values.csv'

enzyme_names = ['topoI', 'RNAP']

title = 'RNAP Tracking Optimization for dt='+str(dt) + ' and '+str(int(percentage_threshold*100))+'% of best cases'

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16
ms=5

# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we load
fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True, sharex=True)

ax = axs

df = pd.read_csv(path + enzyme_names[0]+loss_file)
df = df.sort_values(by='loss', ascending=False)#, inplace=True)
n = len(df['loss'])
nconsidered = int(n*percentage_threshold)
err_threshold = df['loss'].iloc[-nconsidered]
print('Number of tests', n)
print('Considered', nconsidered)
print('For ', percentage_threshold*100, '%')
# Filter according error
filtered_df = df[df['loss'] <= err_threshold]

loss = df['loss'].to_numpy()
floss = filtered_df['loss'].to_numpy()
ax.set_title(title)
ax.plot(df['test'], df['loss'], 'o', ms=ms, color='blue', label='all')
ax.plot(filtered_df['test'], filtered_df['loss'], 'o', ms=ms, color='red', label='best')
#ax.plot(df['loss'], 'o', ms=ms, color='blue')
#ax.plot(loss, 'o', ms=ms, color='blue')
#ax.plot(floss, 'o', ms=ms, color='green')

ax.grid(True)
ax.set_xlabel('test')
ax.set_ylabel('loss')
ax.set_yscale('log')

# Let's calculate averages and save the info:
#-----------------------------------------------------------------------------------------------------------------------
# We will store the dataframes in these lists so we can combine them at the end to create a table.
avg_list = []
std_list = []
topo_clist = ['RNAP_dist', 'fold_change']
RNAP_clist = ['gamma', 'k_on', 'k_off', 'k_open', 'k_closed', 'k_ini']
my_clist = [topo_clist, RNAP_clist]

for i, name in enumerate(enzyme_names):
    df = pd.read_csv(path + name + loss_file)
    df = df.sort_values(by='loss', ascending=False)  # , inplace=True)
    n = len(df['loss'])
    nconsidered = int(n * percentage_threshold)
    err_threshold = df['loss'].iloc[-nconsidered]

    # Filter according error
    filtered_df = df[df['loss'] <= err_threshold]

    # Drop the loss and test columns
    filtered_df = filtered_df.drop(columns=['loss', 'test'])

    # Calculate averages and standard deviations
    df_avg = filtered_df.mean(axis=0).to_frame().T.rename(index={0:'avg'})
    df_std = filtered_df.std(axis=0).to_frame().T.rename(index={0:'std'})

    # Save averages so we can load them for running
    df_avg.to_csv(param_out+ name+ '_dt' + str(dt) + '.csv', index=False, sep=',')

    # Remove columns we don't want for the table
    df_avg = df_avg.loc[:, df_avg.columns.isin(my_clist[i])]
    df_std = df_std.loc[:, df_std.columns.isin(my_clist[i])]

    # Add them to a list
    avg_list.append(df_avg)
    std_list.append(df_std)

# Join them for the table
df_avg = pd.concat([avg_list[0], avg_list[1]], axis=1)
df_std = pd.concat([std_list[0], std_list[1]], axis=1)
new_df = pd.concat([df_avg, df_std], axis=0)
new_df.to_csv(param_out+'table_dt'+str(dt)+'.csv', index=False, sep=',')
#new_df.to_csv('table_dt'+str(dt)+'.csv', index=False, sep=',')

plt.show()
