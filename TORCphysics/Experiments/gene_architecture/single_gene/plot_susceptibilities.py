import pandas as pd
import matplotlib.pyplot as plt

# Description
# ----------------------------------------------------------------------------------------------------------------------
# Having produced results for the biophysical model and collected the relevant experimental data, here we
# collect the corresponding results and output them so we can quickly load them and process them.
# Additionally, we try to plot a graph similar to the one showed in their paper.

# Inputs
# ----------------------------------------------------------------------------------------------------------------------
# TORCphysics data
weak_torc = 'torcphys/production_rates-weak.csv'
medium_torc = 'torcphys/production_rates-medium.csv'
strong_torc = 'torcphys/production_rates-strong.csv'

# Junier experimental data
weak_exp = 'junier_data/weak.csv'
medium_exp = 'junier_data/medium.csv'
strong_exp = 'junier_data/strong.csv'

# Path to biophysical model results
model_weak = 'junier_data/model_weak-norm.csv'
model_medium = 'junier_data/model_medium-norm.csv'
model_strong = 'junier_data/model_strong-norm.csv'

# Read each csv file
weak_torc = pd.read_csv(weak_torc)
medium_torc = pd.read_csv(medium_torc)
strong_torc = pd.read_csv(strong_torc)

weak_exp = pd.read_csv(weak_exp)
medium_exp = pd.read_csv(medium_exp)
strong_exp = pd.read_csv(strong_exp)

model_weak = pd.read_csv(model_weak)
model_medium = pd.read_csv(model_medium)
model_strong = pd.read_csv(model_strong)

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

# Plot
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True, sharex=True)

# Plot experimental curves
axs[0].plot(weak_exp['distance'], weak_exp['Signal'], 's', ms=ms, color=exp_color)
axs[1].plot(medium_exp['distance'], medium_exp['Signal'], 's', ms=ms, color=exp_color)
axs[2].plot(strong_exp['distance'], strong_exp['Signal'], 's', ms=ms, color=exp_color, label='experimental')

# Plot biophysical model curves
axs[0].plot(model_weak['distance'], model_weak['mean'], '--', color=model_color)
axs[1].plot(model_medium['distance'], model_medium['mean'], '--', color=model_color)
axs[2].plot(model_strong['distance'], model_strong['mean'], '--', color=model_color, label='model')

# Plot torcphys curves
axs[0].plot(weak_torc['distance'], weak_torc['prod_rate']/weak_torc['prod_rate'][4], '-', color=torc_color)
axs[1].plot(medium_torc['distance'], medium_torc['prod_rate']/medium_torc['prod_rate'][4], '-', color=torc_color)
axs[2].plot(strong_torc['distance'], strong_torc['prod_rate']/strong_torc['prod_rate'][4], '-', color=torc_color, label='TORC')


for i in range(3):
    ax = axs[i]
    ax.set_xscale('log')
    ax.grid(True)
axs[0].set_ylabel('Weak promoter', fontsize=xlabel_size)
axs[1].set_ylabel('Medium promoter', fontsize=xlabel_size)
axs[2].set_ylabel('Strong promoter', fontsize=xlabel_size)
axs[2].set_xlabel('Distance (bp)', fontsize=xlabel_size)
axs[2].legend(loc='best', fontsize=font_size)

plt.savefig('susceptibility.png')
#plt.savefig('junier.pdf')
plt.show()

