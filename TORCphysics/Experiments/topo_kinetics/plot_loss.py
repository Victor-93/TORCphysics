import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Description
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot the losses
# TODO: What if we plot losses, get avgs/stds using the threshold that at least includes 10% of the
#  best tests, then let's do the optimization plot but for those averages, and see what it gives

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
percentage_threshold = .10
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0 #0.25
initial_time = 0
final_time = 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'loss-'+str(dt)
#file_out = 'calibration_small_dt1'

loss_file = 'recognition-linear/calibration_dt'+str(dt)+'_values.csv'

# TODO: AQUIMEQUEDE
title = 'Topoisomerase optimization for dt='+str(dt) + ' and '+str(int(percentage_threshold*100))+'% of best cases'

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

df = pd.read_csv(loss_file)
df = df.sort_values(by='loss', ascending=False)#, inplace=True)
n = len(df['loss'])
nconsidered = int(n*percentage_threshold)
err_threshold = df['loss'].iloc[-nconsidered]

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

# Calculate averages and standard deviations
df_avg = filtered_df.mean(axis=0).to_frame().T.rename(index={0:'avg'})
df_std = filtered_df.std(axis=0).to_frame().T.rename(index={0:'std'})

# Join them for the table
new_df = pd.concat([df_avg, df_std], axis=0)
new_df.to_csv('table_dt'+str(dt)+'.csv', index=False, sep=',')

# And let's separate them into topoI and gyrase so we can load them.

# topoI
topos_names = ['topoI', 'gyrase']
for name in topos_names:
    topo_dict = {}
    topo_dict['k_on'] = df_avg['k_on_'+name]
    topo_dict['k_off'] = df_avg['k_off_' + name]
    topo_dict['k_cat'] = df_avg['k_cat_'+name]
    topo_dict['k_width'] = df_avg['width_'+name]
    topo_dict['k_threshold'] = df_avg['threshold_'+name]
    if name == 'gyrase':
        topo_dict['sigma0'] = df_avg['sigma0_'+name]

    topo_df = pd.DataFrame.from_dict(topo_dict)
    topo_df.to_csv(name+'_rec_avg_dt'+str(dt)+'.csv', index=False, sep=',')
df_avg.to_csv('avg_dt'+str(dt)+'.csv', index=False, sep=',')

plt.show()
