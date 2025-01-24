import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns


# Description
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot the losses for each case and a distribution of losses at the end

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
percentage_threshold = 0.05

# This list of dicts will include all the input information
info_list = [
    # genearch_V2
    {**{pm: f"genearch_V2/sus_genearch_V2-01-{pm}_dt1.0-trials.pkl" for pm in promoter_cases}, 'label':'V2'},

    # genearch_V2_small-distance
    {**{pm: f"genearch_V2_small-distance/sus_genearch_V2-01-small-dist-{pm}_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V2_sd'},

    # genearch_V0 - gamma 0.01
#    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.01_dt1.0-trials.pkl" for pm in promoter_cases},
#     'label': 'V0_g0.01'},
    # genearch_V0 - gamma 0.05
    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.05_dt1.0-trials.pkl" for pm in promoter_cases},
     'label': 'V0_g0.05'},
    # genearch_V0 - gamma 0.1
#    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.1_dt1.0-trials.pkl" for pm in promoter_cases},
#     'label': 'V0_g0.1'},
    # genearch_V0 - gamma 0.2
#    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.2_dt1.0-trials.pkl" for pm in promoter_cases},
#     'label': 'V0_g0.2'},
    # genearch_V0 - gamma 0.5
#    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.5_dt1.0-trials.pkl" for pm in promoter_cases},
#     'label': 'V0_g0.5'},

    # genearch_V1 - gamma 0.01
#    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.01_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V1_g0.01'},
    # genearch_V1 - gamma 0.05
#    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.05_dt1.0-trials.pkl" for pm in promoter_cases},'label': 'V1_g0.05'},
    # genearch_V1 - gamma 0.1
    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.1_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V1_g0.1'}#,
    # genearch_V1 - gamma 0.2
#    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.2_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V1_g0.2'},
    # genearch_V1 - gamma 0.5
#    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.5_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V1_g0.5'}
]

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 6
height = 3.75
lw = 4
font_size = 12
xlabel_size = 14
title_size = 16
ms=6

colors = ['green', 'blue', 'red']


for case in info_list:
    for pm in promoter_cases:
        with open(case[pm], 'rb') as file:
            data = pickle.load(file)


# Process and plot
#-----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(len(info_list)+1, 3, figsize=(3*width, (len(info_list)+1)*height),
tight_layout=True)#, sharex=True)

# These dataframes will include the distributions of values within the thresholds
loss_dict = {'weak': pd.DataFrame(), 'medium': pd.DataFrame(), 'strong': pd.DataFrame()}

# Let's first plot the all losses
#-----------------------------------------------------------------------------------------------------------------------
for i, case in enumerate(info_list):

    # Go through promoter
    # ---------------------------------------------------
    for j, pm in enumerate(promoter_cases):

        # Load
        # ---------------------------------------------------
        with open(case[pm], 'rb') as file:
            data = pickle.load(file).results

        # Process and plot losses
        # ---------------------------------------------------
        ax = axs[i,j]
        title = case['label'] + ' ' + pm

        # Adjust number of considered cases so all dists have the same number of values
        n = len(data)
        nconsidered = int(n * percentage_threshold)

        # Filter data to include only entries with 'status' == 'ok'
        results = [item for item in data if item['status'] == 'ok']

        system_loss_df = pd.DataFrame([t['loss'] for t in results])

        loss_df = pd.DataFrame({'loss': [t['loss'] for t in results]})

        # Assuming loss_df is a single-column DataFrame
        system_loss_df['loss'] = loss_df.squeeze()  # Convert loss_df to Series if needed

        system_loss_df = system_loss_df.sort_values(by='loss', ascending=False)  # , inplace=True)

        err_threshold = system_loss_df['loss'].iloc[-nconsidered]
        print('Number of tests', n)
        print('Considered', nconsidered)
        print('For ', percentage_threshold * 100, '%')
        # Filter according error
        filtered_df = system_loss_df[system_loss_df['loss'] <= err_threshold]


        loss = system_loss_df['loss'].to_numpy()
        x = np.arange(1, len(loss) + 1, 1)
        ax.plot(x, loss, 'o', ms=ms, color='blue')
        ax.plot(x[n - nconsidered:], loss[n - nconsidered:], 'o', ms=ms, color='red')

        ax.grid(True)
        ax.set_title(title)
        ax.set_xlabel('test')
        ax.set_ylabel('loss')
        ax.set_yscale('log')

        # Add values to loss_dict for plotting distributions
        loss_dict[pm][case['label']] = filtered_df['loss'].reset_index(drop=True)

        # Find the dictionary with the smallest 'loss'
        # data = min(data_list, key=lambda x: x['loss'])
        # x = data['distance']

# Let's now plot the distributions
# -----------------------------------------------------------------------------------------------------------------------
# Go through promoter
# ---------------------------------------------------
for j, pm in enumerate(promoter_cases):
    ax = axs[len(info_list), j]
    title = 'Distributions for ' +str(100*percentage_threshold) + '% of best - ' + pm

    sns.violinplot(data=loss_dict[pm], ax=ax, inner="quart")  # , cut=0, color=colors[i])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel('Loss')
    ax.grid(True)
    ax.set_title(title)

#plt.savefig('V0-V1-V2_losses.png')
plt.show()

