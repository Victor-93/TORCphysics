import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import seaborn as sns


# Description
#-----------------------------------------------------------------------------------------------------------------------
#  Let's plot all general information  as a distribution, where each row is a variable and columns are weak,
#       medium and strong promoters. For each plot, we put together all cases V1/V2. Variables include loss, prod_rate,
#       local/global supercoiling, avg number of bound RNAPs, topos, and gyrases
#       This script has similar structura than analyse_trials-loss.py

# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
variables = ['loss', 'prod_rate', 'local_superhelical', 'global_superhelical']
percentage_threshold = 0.05

# This list of dicts will include all the input information
info_list = [
    # genearch_V2
    {**{pm: f"genearch_V2/sus_genearch_V2-01-{pm}_dt1.0-trials.pkl" for pm in promoter_cases}, 'label':'V2'},

    # genearch_V0 - gamma 0.01
    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.01_dt1.0-trials.pkl" for pm in promoter_cases},
     'label': 'V0_g0.01'},
    # genearch_V0 - gamma 0.05
    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.05_dt1.0-trials.pkl" for pm in promoter_cases},
     'label': 'V0_g0.05'},
    # genearch_V0 - gamma 0.1
    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.1_dt1.0-trials.pkl" for pm in promoter_cases},
     'label': 'V0_g0.1'},
    # genearch_V0 - gamma 0.2
    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.2_dt1.0-trials.pkl" for pm in promoter_cases},
     'label': 'V0_g0.2'},
    # genearch_V0 - gamma 0.5
    {**{pm: f"genearch_V0/sus_genearch_V0-01-{pm}_gamma0.5_dt1.0-trials.pkl" for pm in promoter_cases},
     'label': 'V0_g0.5'},

    # genearch_V1 - gamma 0.01
    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.01_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V1_g0.01'},
    # genearch_V1 - gamma 0.05
    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.05_dt1.0-trials.pkl" for pm in promoter_cases},'label': 'V1_g0.05'},
    # genearch_V1 - gamma 0.1
    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.1_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V1_g0.1'},
    # genearch_V1 - gamma 0.2
    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.2_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V1_g0.2'},
    # genearch_V1 - gamma 0.5
    {**{pm: f"genearch_V1/sus_genearch_V1-01-{pm}_gamma0.5_dt1.0-trials.pkl" for pm in promoter_cases}, 'label': 'V1_g0.5'}
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
#fig, axs = plt.subplots(len(variables), 3, figsize=(3*width, len(variables)*height),
fig, axs = plt.subplots(7, 3, figsize=(3*width, 7*height),
tight_layout=True)#, sharex=True)

# Let's plot as we load

# Go through promoter
# ---------------------------------------------------
for j, pm in enumerate(promoter_cases):

#    ax = axs[v, j]
 #   ylabel = 'loss'
 #   title = ylabel + ' ' + pm

    #var_df = pd.DataFrame() # We will collect data here and plot it
    loss_df = pd.DataFrame()
    prod_rate_df = pd.DataFrame()
    local_superhelical_df = pd.DataFrame()
    global_superhelical_df = pd.DataFrame()
    nRNAPs_df = pd.DataFrame()
    ntopoIs_df = pd.DataFrame()
    ngyrases_df = pd.DataFrame()

    # Go through case
    # ---------------------------------------------------
    for i, case in enumerate(info_list):

        # Load
        # ---------------------------------------------------
        with open(case[pm], 'rb') as file:
            data = pickle.load(file).results

        # ---------------------------------------------------------------------------------
        # Collect variables
        # ---------------------------------------------------------------------------------
        # Adjust number of considered cases so all dists have the same number of values
        n = len(data)
        nconsidered = int(n * percentage_threshold)

        # Filter data to include only entries with 'status' == 'ok'
        results = [item for item in data if item['status'] == 'ok']

        # Collect variables
        loss = [t['loss'] for t in results]
        prod_rate = [np.mean(t['prod_rate']) for t in results]
        local_superhelical = [np.mean(t['local_superhelical']) for t in results]
        global_superhelical = [np.mean(t['global_superhelical']) for t in results]
        nRNAPs = [np.mean(t['nenzymes']['RNAP']) for t in results]
        ntopoIs = [np.mean(t['nenzymes']['topoI']) for t in results]
        ngyrases = [np.mean(t['nenzymes']['gyrase']) for t in results]

        # And put them in df
        case_df = pd.DataFrame({
            "loss": loss,
            "prod_rate": prod_rate,
            "local_superhelical": local_superhelical,
            "global_superhelical": global_superhelical,
            "nRNAPs": nRNAPs,
            "ntopoIs": ntopoIs,
            "ngyrases": ngyrases
        })

        # ---------------------------------------------------------------------------------
        # Filter according loss
        # ---------------------------------------------------------------------------------
        case_df = case_df.sort_values(by='loss', ascending=False)  # , inplace=True)
        err_threshold = case_df['loss'].iloc[-nconsidered]
        print('Number of tests', n)
        print('Considered', nconsidered)
        print('For ', percentage_threshold * 100, '%')
        # Filter according error
        filtered_df = case_df[case_df['loss'] <= err_threshold]

        # ---------------------------------------------------------------------------------
        # Add values to dataframes
        # ---------------------------------------------------------------------------------
        loss_df[case['label']] = filtered_df['loss'].reset_index(drop=True)
        prod_rate_df[case['label']] = filtered_df['prod_rate'].reset_index(drop=True)
        local_superhelical_df[case['label']] = filtered_df['local_superhelical'].reset_index(drop=True)
        global_superhelical_df[case['label']] = filtered_df['global_superhelical'].reset_index(drop=True)
        nRNAPs_df[case['label']] = filtered_df['nRNAPs'].reset_index(drop=True)
        ntopoIs_df[case['label']] = filtered_df['ntopoIs'].reset_index(drop=True)
        ngyrases_df[case['label']] = filtered_df['ngyrases'].reset_index(drop=True)


    # ---------------------------------------------------------------------------------
    # Plot variables!
    # ---------------------------------------------------------------------------------
    # Loss
    ax = axs[0, j]
    ylabel = 'loss'
    title = ylabel + ' for t ' +str(percentage_threshold) + ' ' + pm

    sns.violinplot(data=loss_df, ax=ax, inner="quart")  # , cut=0, color=colors[i])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.set_title(title)
    ax.set_yscale('log')

    # Production rate
    ax = axs[1, j]
    ylabel = 'prod_rate'
    title = ylabel + ' for t ' +str(percentage_threshold) + ' ' + pm

    sns.violinplot(data=prod_rate_df, ax=ax, inner="quart")  # , cut=0, color=colors[i])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.set_title(title)

    # Local superhelical
    ax = axs[2, j]
    ylabel = 'Local superhelical'
    title = ylabel + ' for t ' + str(percentage_threshold) + ' ' + pm

    sns.violinplot(data=local_superhelical_df, ax=ax, inner="quart")  # , cut=0, color=colors[i])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.set_title(title)

    # Global superhelical
    ax = axs[3, j]
    ylabel = 'Global superhelical'
    title = ylabel + ' for t ' + str(percentage_threshold) + ' ' + pm

    sns.violinplot(data=global_superhelical_df, ax=ax, inner="quart")  # , cut=0, color=colors[i])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.set_title(title)

    # nRNAPs
    ax = axs[4, j]
    ylabel = 'RNAPs'
    title = ylabel + ' for t ' + str(percentage_threshold) + ' ' + pm

    sns.violinplot(data=nRNAPs_df, ax=ax, inner="quart")  # , cut=0, color=colors[i])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.set_title(title)

    # nTopoIs
    ax = axs[5, j]
    ylabel = 'TopoIs'
    title = ylabel + ' for t ' + str(percentage_threshold) + ' ' + pm

    sns.violinplot(data=ntopoIs_df, ax=ax, inner="quart")  # , cut=0, color=colors[i])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.set_title(title)

    # ngyrases
    ax = axs[6, j]
    ylabel = 'Gyrases'
    title = ylabel + ' for t ' + str(percentage_threshold) + ' ' + pm

    sns.violinplot(data=ngyrases_df, ax=ax, inner="quart")  # , cut=0, color=colors[i])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.set_title(title)

plt.show()
