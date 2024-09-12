import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle


# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
dt=0.5

experimental_files = []
calibration_files = []
for pcase in promoter_cases:
    experimental_files.append('../../junier_data/inferred-rate_' + pcase + '.csv')
    #calibration_files.append('Stages-'+pcase+'_dt1.pkl')
    calibration_files.append('Stages-'+pcase+'_dt'+str(dt)+'.pkl')

# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 8
height = 4
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

# line styles
model_ls = '-'
exp_ls = '--'
titles = ['weak', 'medium', 'strong']

colors = ['green', 'blue', 'red']

# Processing functions
#-----------------------------------------------------------------------------------------------------------------------
# Calculate rates as a function of distance
def get_prod_rates(results_list):
    x = [d['distance'] for d in results_list]
    # y = np.zeros_like(x)
    # ys = np.zeros_like(x)
    y = []
    ys = []
    # List with just results (not distance and not dict anymore)
    results_list2 = [d['result']['prod_rate'] for d in results_list]
    for j, rates_array in enumerate(results_list2): #case_results is the rate
        #rates = [d[0] for d in case_results]

        # Convert to a NumPy array
        #rates_array = np.array(case_re)

        # Calculate mean and standard deviation
        mean = np.mean(rates_array)
        std = np.std(rates_array)
        #y[j] = mean
        #ys[j] = std
        y.append(mean)
        ys.append(std)

    rates = np.array([x, y, ys])
    return rates

# Load
#-----------------------------------------------------------------------------------------------------------------------
pickle_data = []
for pickle_file in calibration_files:
    with open(pickle_file, 'rb') as file:
        data = pickle.load(file)
        pickle_data.append(data)

# Calculate rates
rates = []
for i, data in enumerate(pickle_data):
    x=i
    rates.append(get_prod_rates(data[0]['data']))

# Plot
#-----------------------------------------------------------------------------------------------------------------------
# Let's plot as we do the process
fig, axs = plt.subplots(4, figsize=(width, 4*height), tight_layout=True, sharex=True)

for i, rate_array in enumerate(rates):

    axs[i].set_title(titles[i])

    # Prepare calibration array and plot
    x = rate_array[0]
    y = rate_array[1]
    ys = rate_array[2]

    axs[i].plot(x, y, model_ls, lw=lw, color=colors[i])
    axs[3].plot(x, y, model_ls, lw=lw, color=colors[i])     # Last one

    # Load experimental and plot
    exp = pd.read_csv(experimental_files[i]) # read
    x = exp['distance']
    y = exp['Signal']
    ys = exp['Error']
    axs[i].plot(x, y, exp_ls, lw=lw, color=colors[i])
    axs[3].plot(x, y, exp_ls, lw=lw, color=colors[i])     # Last one

    axs[i].set_xlabel('Distance')
    axs[i].set_ylabel('expression rate')
    axs[i].grid(True)
    axs[i].set_xscale('log')


# Last one
axs[3].set_title('all together')
axs[3].set_xlabel('Distance')
axs[3].set_ylabel('expression rate')
axs[3].grid(True)
axs[3].set_xscale('log')

plt.show()



