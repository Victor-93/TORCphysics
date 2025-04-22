import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Description
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot the losses

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
percentage_threshold = 0.05#005#.10
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0 #0.25
initial_time = 0
final_time = 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'topoI-tracking-loss-'+str(dt)
#file_out = 'TORCMeeting-tracking-loss-'+str(dt)
path = 'track-StagesStall/avg02/'
loss_file = '-calibration_avg-RNAPTracking_nsets_p2_small_dt'+str(dt)+'-values.csv'

enzyme_names = ['topoI', 'RNAP']

title = 'RNAP Tracking Optimization for dt='+str(dt) + ' and '+str(int(percentage_threshold*100))+'% of best cases'
title = 'TopoI-RNAP tracking calibration'

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

df = pd.read_csv(path+enzyme_names[0]+loss_file)
df = df.sort_values(by='loss', ascending=False)#, inplace=True)
n = len(df['loss'])
nconsidered = int(n*percentage_threshold)
err_threshold = df['loss'].iloc[-nconsidered]

# Filter according error
filtered_df = df[df['loss'] <= err_threshold]

loss = df['loss'].to_numpy()
floss = filtered_df['loss'].to_numpy()
ax.set_title(title)

# Create a histogram
minv = min(loss)
maxv = np.mean(loss) + 1*np.std(loss)
#maxv = max(loss)*.0075
maxv = .175#max(loss)*.2
bins = np.linspace(minv, maxv, 150)  # Define bins
hist, bin_edges = np.histogram(loss, bins=bins)

# Plot the full histogram
ax.hist(loss, bins=bins, color='gray', alpha=0.6, label='data')

selected_loss =floss
#selected_loss = [floss[-1]]
# Highlight bins corresponding to floss
i=1
for value in selected_loss:
    # Find the bin index for the current value
    bin_index = np.digitize(value, bin_edges) - 1
    # Plot the specific bin
    plt.bar(
        bin_edges[bin_index],  # Bin start
        hist[bin_index],  # Bin height
        width=bin_edges[1] - bin_edges[0],  # Bin width
        color='red',  # Highlight color
        alpha=0.8,
#        edgecolor='black',
        label='filtered' if i == len(selected_loss) else ""
    )
    i = i + 1
ax.set_title(title, fontsize=title_size)
ax.set_xlabel('loss $\epsilon$', fontsize=xlabel_size)
ax.set_ylabel(r'$n$ tests', fontsize=xlabel_size)
ax.legend(loc='best')

#ax.set_xscale('log')
ax.grid(True)

plt.savefig(file_out+'.pdf', dpi=300)
#plt.savefig(file_out+'.png', dpi=300)
plt.show()
