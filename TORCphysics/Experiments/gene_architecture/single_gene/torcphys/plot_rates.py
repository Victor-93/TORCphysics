import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle


# Inputs
#-----------------------------------------------------------------------------------------------------------------------
promoter_cases = ['weak', 'medium', 'strong']
dt=1.0#0.5

k_weak=0.03
# k_weak=0.0334

model_code = 'GB-Stages-'

experimental_files = []
calibration_files = []
for pcase in promoter_cases:
    # experimental_files.append('../../junier_data/inferred-rate_' + pcase + '.csv')
    experimental_files.append('../junier_data/inferred-rate_kw'+str(k_weak)+'_' + pcase + '.csv')
    #calibration_files.append('Stages-'+pcase+'_dt1.pkl')
    #calibration_files.append('Stages-'+pcase+'_dt'+str(dt)+'.pkl')
    #calibration_files.append('k_ini-Stages-'+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.pkl')
    # calibration_files.append('GB-Stages-'+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.pkl')

    # calibration_files.append('calibrate_inferred-rates/'+model_code+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.pkl')
    calibration_files.append('calibrate_inferred-rates/reproduce-'+model_code+pcase+'-kw'+str(k_weak)+'_dt'+str(dt)+'.pkl')


# Plotting params
#-----------------------------------------------------------------------------------------------------------------------
width = 7
height = 3
lw = 3
font_size = 12
xlabel_size = 14
title_size = 16

# line styles
model_ls = '-o'
exp_ls = '--o'
titles = ['Weak promoter', 'Medium promoter', 'Strong promoter']

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
        mean = rates_array[0]
        std = rates_array[1]
        #mean = np.mean(rates_array)
        #std = np.std(rates_array)
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

    axs[i].set_title(titles[i], color=colors[i])

    # Prepare calibration array and plot
    x = rate_array[0]
    y = rate_array[1]
    ys = rate_array[2]#/np.sqrt(n_sims)

    if i == 0:
        axs[i].plot(x, y, model_ls, lw=lw, color=colors[i], label='TORCPhysics')
    else:
        axs[i].plot(x, y, model_ls, lw=lw, color=colors[i])

    axs[i].fill_between(x, y-ys, y+ys, lw=lw, color=colors[i], alpha=0.2)
    axs[3].plot(x, y, model_ls, lw=lw, color=colors[i])     # Last one
    axs[3].fill_between(x, y-ys, y+ys, lw=lw, color=colors[i], alpha=0.1)

    # Load experimental and plot
    exp = pd.read_csv(experimental_files[i]) # read
    x = exp['distance']
    y = exp['Signal']
    ys = exp['Error']
    if i == 0:
        # axs[i].plot(x, y, exp_ls, lw=lw, color=colors[i], label='experiment')
        axs[i].errorbar(x, y, yerr=ys, fmt=exp_ls, capsize=5,color=colors[i], label='Experiment')
    else:
        #axs[i].plot(x, y, exp_ls, lw=lw, color=colors[i])
        axs[i].errorbar(x, y, yerr=ys, fmt=exp_ls, capsize=5,color=colors[i])

    #axs[3].plot(x, y, exp_ls, lw=lw, color=colors[i])     # Last one
    axs[3].errorbar(x, y, yerr=ys, fmt=exp_ls, capsize=5, color=colors[i])

    # Error quantification - meansquared error
    # error = np.sum((exp['Signal'] - rate_array[1]) ** 2)#/len(rate_array[1])
    error = np.sum((exp['Signal'] - rate_array[1]) ** 2)/len(rate_array[1])
    formatted_sum_err = "{:.3g}".format(error)  # Formatting to three significant figures
    er_label = r'$\epsilon={}$'.format(formatted_sum_err)

    # Add the text box to the plot
    props = dict(boxstyle='round', facecolor='silver', alpha=0.5)
    axs[i].text(0.7, 0.1, er_label, transform=axs[i].transAxes, fontsize=14,
            verticalalignment='top', bbox=props)

    # axs[i].set_xlabel('Upstream distance (bp)')
    axs[i].set_ylabel(r'Expression rate $(s^{-1})$')
    axs[i].grid(True)
    axs[i].set_xlim([80, 6000])
    axs[i].set_xscale('log')

axs[3].set_xlim([80, 6000])

# Last one
axs[3].set_title('All together')
axs[3].set_xlabel('Upstream distance (bp)')
axs[3].set_ylabel(r'Expression rate $(s^{-1})$')
axs[3].grid(True)
axs[3].set_xscale('log')
axs[0].legend(loc='best')

#plt.savefig(model_code+'kw'+str(k_weak)+'_dt'+str(dt)+'.png')
#plt.savefig(model_code+'kw'+str(k_weak)+'_dt'+str(dt)+'.pdf')

plt.show()



