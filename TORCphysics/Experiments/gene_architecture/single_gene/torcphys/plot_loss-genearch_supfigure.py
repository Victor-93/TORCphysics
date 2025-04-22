import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Description
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot the losses

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
promoter_labels = ['Weak', 'Medium', 'Strong']
percentage_threshold = 0.004 # Because we choose the best one here, not a set
# Units:
# concentrations (nM), K_M (nM), velocities (nM/s), time (s)
dt = 1.0 #0.25
initial_time = 0
final_time = 500
time = np.arange(initial_time, final_time + dt, dt)
frames = len(time)
file_out = 'genearch-loss-'+str(dt)

# This list of dicts will include all the input information
loss_files = [
    # genearch_V0 - gamma 0.05
    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.05_dt1.0-values.csv" for pm in promoter_cases},
     'label': 'V0'},

    # genearch_V1 - gamma 0.1
    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.1_dt1.0-values.csv" for pm in promoter_cases}, 'label': 'V1'},

    # genearch_V2
    {**{pm: f"genearch_V2/sus_genearch_V2-01-{pm}_dt1.0-values.csv" for pm in promoter_cases}, 'label': 'V2'}
]

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 6#8
height = 3#4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16
ms=5

colors = ['green', 'blue', 'red', 'purple', 'orange']

# Plot
#-----------------------------------------------------------------------------------------------------------------------
ncols = len(loss_files)
nrows = len(promoter_cases)
fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*width, nrows*height),
tight_layout=True)

ranges = [0.5, 0.5, 0.5]
# ---------------------------------------------------
s= -1
for i, case in enumerate(loss_files):

    # Go through promoter
    for j, pm in enumerate(promoter_cases):

        s=s+1 # This is for the outside label

        ax = axs[j, i]
        loss_file = case[pm]

        figtitle = case['label'] + ': ' + promoter_labels[j] + ' Promoter'
        ax.set_title(figtitle, fontsize=title_size)  # color=colors[i],

        df = pd.read_csv(loss_file)
        df = df.sort_values(by='loss', ascending=False)#, inplace=True)
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df = df.dropna()
        n = len(df['loss'])
        nconsidered = int(n*percentage_threshold)
        err_threshold = df['loss'].iloc[-nconsidered]

        # Filter according error
        filtered_df = df[df['loss'] <= err_threshold]

        loss = df['loss'].to_numpy()
        floss = filtered_df['loss'].to_numpy()

        # Create a histogram
        minv = min(loss)
        #maxv = np.mean(loss) + 1*np.std(loss)
        #maxv = 0.5#1.#max(loss)*.5
        #maxv = max(loss)*.1
        maxv = ranges[i]
        bins = np.linspace(minv, maxv, 100)  # Define bins
        hist, bin_edges = np.histogram(loss, bins=bins)

        # Plot the full histogram
        ax.hist(loss, bins=bins, color='gray', alpha=0.6, label='data')

        selected_loss =floss
        #selected_loss = [floss[-1]]
        # Highlight bins corresponding to floss
        k=1
        for value in selected_loss:
            # Find the bin index for the current value
            bin_index = np.digitize(value, bin_edges) - 1
            # Plot the specific bin
            ax.bar(
                bin_edges[bin_index],  # Bin start
                hist[bin_index],  # Bin height
                width=bin_edges[1] - bin_edges[0],  # Bin width
                color='red',  # Highlight color
                alpha=0.8,
        #        edgecolor='black',
                label='filtered' if k == len(selected_loss) else ""
            )
            k=k+1
        ax.set_xlabel('loss $\epsilon$', fontsize=xlabel_size)
        ax.set_ylabel(r'$n$ tests', fontsize=xlabel_size)
        ax.grid(True)
axs[0,0].legend(loc='best')

plt.savefig(file_out+'.pdf', dpi=300)
plt.show()
