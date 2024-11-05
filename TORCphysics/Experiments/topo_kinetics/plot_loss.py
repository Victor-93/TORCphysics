import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Description
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot the losses

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
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
titles = ['Weak dt=1.0', 'Medium dt=1.0', 'Strong dt=1.0']

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
#df = df.sort_values(by='loss', ascending=False)#, inplace=True)
#df.sort_values(by=['loss'], inplace=True)
loss = df['loss'].to_numpy()
#ax.set_title(titles[i])
#ax.plot(df['test'], df['loss'], 'o', ms=ms, color='blue')
#ax.plot(df['loss'], 'o', ms=ms, color='blue')
ax.plot(loss, 'o', ms=ms, color='blue')

# Let's get the 5 smallest values:
#df.sort_values(by=['loss'], inplace=True)
#mdf = df.nsmallest(5, 'loss')
#ax.plot(mdf['test'], mdf['loss'], 'o', ms=ms*1.5, color='red')
#ax.plot(df['test'].iloc[0], mdf['loss'].iloc[0], 'o', ms=ms*2, color='green') # And the best

ax.grid(True)
ax.set_xlabel('test')
ax.set_ylabel('loss')
ax.set_yscale('log')

# Let's print some info:
#print(titles[i])
#print(mdf)


plt.show()
