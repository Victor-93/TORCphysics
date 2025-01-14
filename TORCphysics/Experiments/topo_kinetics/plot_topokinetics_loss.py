import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Description
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot the losses
# TODO: Keep working on it

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

# Create a histogram
minv = min(loss)
maxv = np.mean(loss) + 1*np.std(loss)
maxv = .5#max(loss)*.2
bins = np.linspace(minv, maxv, 100)  # Define bins
hist, bin_edges = np.histogram(loss, bins=bins)

# Plot the full histogram
ax.hist(loss, bins=bins, color='gray', alpha=0.6, label='Loss')

# Highlight bins corresponding to floss
for value in floss:
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
        label='Highlighted' if bin_index == np.digitize(floss[0], bin_edges) - 1 else ""
    )
ax.grid(True)
#ax.set_xlabel('test')
#ax.set_ylabel('loss')
ax.set_xscale('log')


plt.show()
