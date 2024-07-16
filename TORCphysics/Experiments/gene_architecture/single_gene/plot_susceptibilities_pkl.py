import pandas as pd
import matplotlib.pyplot as plt
import pickle
import numpy as np
import seaborn as sns

# Description
# ----------------------------------------------------------------------------------------------------------------------
# Similar to the original version. But this time I plot from the pickle file.

# Inputs
# ----------------------------------------------------------------------------------------------------------------------
# TORCphysics data
weak_torc = 'torcphys/production_rates-weak.csv'
medium_torc = 'torcphys/production_rates-medium.csv'
strong_torc = 'torcphys/production_rates-strong.csv'

# Pickle file
weak_torc_pkl = 'torcphys/production_rates-weak.pkl'

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

# Read pickle file
pickle_file = weak_torc_pkl
with open(pickle_file, 'rb') as file:
    output_weak = pickle.load(file)

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


# Calculate rates as a function of distance
def get_prod_rates(results_list):
    x = [d['distance'] for d in results_list]
    # y = np.zeros_like(x)
    # ys = np.zeros_like(x)
    y = []
    ys = []
    # List with just results (not distance and not dict anymore)
    results_list2 = [d['results'] for d in results_list]
    for j, case_results in enumerate(results_list2):
        rates = [d[0] for d in case_results]

        # Convert to a NumPy array
        rates_array = np.array(rates)

        # Calculate mean and standard deviation
        mean = np.mean(rates_array)
        std = np.std(rates_array)
        #y[j] = mean
        #ys[j] = std
        y.append(mean)
        ys.append(std)

    rates = np.array([x, y, ys])
    return rates


# Get's the KDEs of RNAPs
def get_RNAP_kde(results_list):
    RNAP_KDE = []
    # List with just results (not distance and not dict anymore)
    results_list2 = [d['results'] for d in results_list]
    for j, case_results in enumerate(results_list2):
        RNAP_kde_dict = [d[1]['KDE']['RNAP'] for d in case_results]
        kde_y = np.array([d['kde_y'] for d in RNAP_kde_dict])

        mean = np.mean(kde_y, axis=0)
        std = np.std(kde_y, axis=0)
        kde = np.array([mean, std])
        RNAP_KDE.append(kde)

    return RNAP_KDE


# Get's global superhelicity, and returns distribution>
def get_global_superhelicity(results_list):
    data = []

    # List with just results (not distance and not dict anymore)
    results_list2 = [d['results'] for d in results_list]
    for j, case_results in enumerate(results_list2):
        dist = results_list[j]['distance']
        sigmas = [d[1]['global_superhelical'] for d in case_results]
        sigmas = [item for sublist in sigmas for item in sublist]
        sigmas = np.array(sigmas)

        data.append({'distance': dist, 'sigma': sigmas})

        # data.append({'distance': np.repeat(dist, len(sigmas)), 'sigma': sigmas})

    # Flatten the data into a list of dictionaries
    flattened_data = []
    for entry in data:
        distance = entry['distance']
        for sigma_value in entry['sigma']:
            flattened_data.append({'distance': distance, 'sigma': sigma_value})
    return pd.DataFrame(flattened_data)


# Get's local superhelicity, and returns distribution>
def get_local_superhelicity(results_list):
    data = []

    # List with just results (not distance and not dict anymore)
    results_list2 = [d['results'] for d in results_list]
    for j, case_results in enumerate(results_list2):
        dist = results_list[j]['distance']
        sigmas = [d[1]['local_superhelical'] for d in case_results]
        sigmas = [item for sublist in sigmas for item in sublist]
        sigmas = np.array(sigmas)

        data.append({'distance': dist, 'sigma': sigmas})

    # Flatten the data into a list of dictionaries
    flattened_data = []
    for entry in data:
        distance = entry['distance']
        for sigma_value in entry['sigma']:
            flattened_data.append({'distance': distance, 'sigma': sigma_value})
    return pd.DataFrame(flattened_data)

# Calculate rates
prod_rates = get_prod_rates(output_weak)

# Calculate RNAP KDEs
RNAP_kdes = get_RNAP_kde(output_weak)

# Get and sort global sigma data sorted by distances so seaborn can create violin plots.
global_sigma = get_global_superhelicity(output_weak)

local_sigma = get_local_superhelicity(output_weak)

kde_x = np.linspace(0, 100, len(RNAP_kdes[0][1]))

# Plot - 1 - rates, 2.- RNAP location - 3.- violin plots of gloval superhelicity
# ----------------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(4, figsize=(width, 4 * height), tight_layout=True, sharex=False)

fig.suptitle('test promoter')

# rates
# --------------------------------------------------------------
# Plot torcphys curves
#axs.plot(prod_rates[0], prod_rates[1]/prod_rates[1][4], '-o', color=torc_color)
axs[0].plot(prod_rates[0], prod_rates[1], '-o', color=torc_color)
axs[0].fill_between(prod_rates[0], prod_rates[1] - prod_rates[2], prod_rates[1] + prod_rates[2], color=torc_color,
                    alpha=0.5)

# RNAP KDEs
# --------------------------------------------------------------

colors = ['red', 'blue', 'green', 'purple', 'black']
indexes = [0, 4, 8, -3, -1]

ax = axs[1]
for count, index in enumerate(indexes):
    axs[0].plot(prod_rates[0][index], prod_rates[1][index], 'o', ms=10, color=colors[count])
    ax.plot(kde_x, RNAP_kdes[index][0], color=colors[count])
    #ax.fill_between(kde_x,
    #                RNAP_kdes[index][0] + RNAP_kdes[index][1],
    #                RNAP_kdes[index][0] - RNAP_kdes[index][1],
    #                color=colors[count], alpha=0.15
    #                )
#for kde in RNAP_kdes:
#    ax.plot(kde_x, kde[1], '-')

# global and local superhelical
# --------------------------------------------------------------
sns.violinplot(x='distance', y='sigma', data=global_sigma, ax=axs[2])
sns.violinplot(x='distance', y='sigma', data=local_sigma, ax=axs[3], color='skyblue')

axs[0].set_xscale('log')
axs[0].grid(True)
axs[0].set_ylabel('Production rate', fontsize=xlabel_size)
axs[0].set_xlabel('upstream distance (bp)', fontsize=xlabel_size)

axs[1].grid(True)
axs[1].set_ylabel('Density', fontsize=xlabel_size)
axs[1].set_xlabel('gene length (%)', fontsize=xlabel_size)

sigma_a = -.08
sigma_b = 0.04
axs[2].grid(True)
axs[2].set_ylabel('Global superhelical level', fontsize=xlabel_size)
axs[2].set_xlabel('upstream distance (bp)', fontsize=xlabel_size)
axs[2].set_ylim(sigma_a, sigma_b)

axs[3].grid(True)
axs[3].set_ylabel('Local superhelical level', fontsize=xlabel_size)
axs[3].set_xlabel('upstream distance (bp)', fontsize=xlabel_size)
axs[3].set_ylim(sigma_a, sigma_b)

#axs.legend(loc='best', fontsize=font_size)

#plt.savefig('susceptibility.png')
#plt.savefig('junier.pdf')
plt.show()
