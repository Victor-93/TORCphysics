import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
loss_files = []
loss_files.append('calibrate_inferred-rates/GB-Stages-weak-kw0.02_dt1.0-values.csv')
loss_files.append('calibrate_inferred-rates/GB-Stages-medium-kw0.02_dt1.0-values.csv')
loss_files.append('calibrate_inferred-rates/GB-Stages-strong-kw0.02_dt1.0-values.csv')

#loss_files.append('calibrate_inferred-rates/Stages-weak_dt1.0-values.csv')
#loss_files.append('calibrate_inferred-rates/Stages-medium_dt1.0-values.csv')
#loss_files.append('calibrate_inferred-rates/Stages-strong_dt1.0-values.csv')

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
fig, axs = plt.subplots(len(loss_files), figsize=(width, len(loss_files)*height), tight_layout=True, sharex=True)

for i, loss_file in enumerate(loss_files):

    if len(loss_files) >1:
        ax = axs[i]
    else:
        ax =axs

    df = pd.read_csv(loss_file)

    ax.set_title(titles[i])
    ax.plot(df['test'], df['loss'], '-o', ms=ms, color='blue')

    # Let's get the 5 smallest values:
    df.sort_values(by=['loss'], inplace=True)
    mdf = df.nsmallest(5, 'loss')
    ax.plot(mdf['test'], mdf['loss'], 'o', ms=ms*1.5, color='red')
    ax.plot(df['test'].iloc[0], mdf['loss'].iloc[0], 'o', ms=ms*2, color='green') # And the best

    ax.grid(True)
    ax.set_xlabel('test')
    ax.set_ylabel('loss')
    ax.set_yscale('log')

    # Let's print some info:
    print(titles[i])
    print(mdf)


plt.show()
